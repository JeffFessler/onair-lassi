%%
% Compare Unitary-DIL and online Unitary-DIL for video denoising
%
% Brian Moore
% brimoor@umich.edu
%

rng(42);

% Data knobs
path = '/Users/Brian/Archive/MATLAB/video_data/road.mat';
idim = [64, nan, 50];               % Raw video dimensions
T    = 5;                           % Streaming size
dt   = 1;                           % Temporal patch stride
pdim = [8, 8, T];                   % Patch sizes
pgap = [2, 2, dt];                  % Patch strides
PROB = 0.05;                        % Missing data probability

% Algorithm knobs
mu       = 0.5;                     % Sparsity regularization
gamma    = 0.9;                     % Forgetting factor (online)
D0       = dctmtx(prod(pdim))';     % Initial dictionary
type     = 'hard';                  % Sparsity regularizer
fixedD   = false;                   % Fixed dictionary?
nIters1  = 15;                      % # Unitary-DIL iters
nIters2i = 15;                      % # iters-init Online Unitary-DIL
nIters2  = 5;                       % # iters/time Online Unitary-DIL
flag     = 1;                       % Print stats?

% Load video
Xtrue = loadVideo(path,idim,T,dt);
[ny, nx, nt] = size(Xtrue);
[p, m] = size(D0);

% Missing data
M = (rand(ny,nx,nt) > PROB);
Y = Xtrue;
Y(~M) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unitary-DIL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patch data
[Mp1, ~, Np1] = patchInds([ny, nx, nt],pdim,pgap);
Z = Y(Mp1);

% Unitary-DIL
opts1 = struct();
opts1.D0     = D0;
opts1.type   = type;
opts1.fixedD = fixedD;
opts1.nIters = nIters1;
opts1.flag   = flag;
[D1, B1, stats1] = unitaryDil(Z,mu,opts1);
cost1 = stats1.cost(:);
repError1 = stats1.repError(:);

% Reconstruction
Xhat1 = patchRecon(D1 * B1,Mp1,Np1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Online Unitary-DIL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
[Mp2, ~, Np2] = patchInds([ny, nx, T],pdim,pgap);
Xhat2 = zeros(ny,nx,nt);
tt = 1:dt:(nt + 1 - T);
ni = numel(tt);
cost2 = nan(nIters2i,ni);
repError2 = nan(nIters2i,ni);

% Patch data
Y1 = Y(:,:,1:T);
Z1 = Y1(Mp2);

% Online Unitary-DIL (t = 1)
opts2 = struct();
opts2.D0     = D0;
opts2.type   = type;
opts2.fixedD = fixedD;
opts2.nIters = nIters2i;
opts2.flag   = flag;
[Dt, Bt, params, stats2] = onlineUnitaryDil(Z1,mu,gamma,opts2);
cost2(:,1) = stats2.cost(:);
repError2(:,1) = stats2.repError(:);

% Initialize reconstruction
Xhat2(:,:,1:T) = patchRecon(Dt * Bt,Mp2,Np2);

% Online Unitary-DIL (t = 2,3,...)
for i = 2:ni
    % Patch data
    t = tt(i);
    Yt = Y(:,:,t:(t + T - 1));
    Zt = Yt(Mp2);
    
    % Online Unitary-DIL
    opts2.D0     = Dt;
    opts2.B0     = Bt;
    opts2.params = params;
    opts2.nIters = nIters2;
    [Dt, Bt, params, stats2] = onlineUnitaryDil(Zt,mu,gamma,opts2);
    cost2(1:nIters2,i) = stats2.cost(:);
    repError2(1:nIters2,i) = stats2.repError(:);
    
    % Update reconstruction
    Xt = patchRecon(Dt * Bt,Mp2,Np2);
    Xhat2 = updateRecon(Xhat2,Xt,gamma,t,dt);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
movie.video = cell2mat({Xtrue,Y;Xhat1,Xhat2});
opts.xlabels = {'Truth','Observations';'Batch','Online'};
PlayMovie(movie,opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cm = linspecer(2);
cfigure();

% Unitary-DIL
subplot(2,2,1);
plot(1:nIters1,cost1,'-o','Color',cm(1,:)); hold on;
xlabel('iteration');
ylabel('cost');
title('unitaryDil');
axis tight; padAxis();

subplot(2,2,3);
plot(1:nIters1,repError1,'-o','Color',cm(1,:)); hold on;
xlabel('iteration');
ylabel('representation error');
title('unitaryDil');
axis tight; padAxis();

% Online Unitary-DIL
subplot(2,2,2);
plot(1:nIters2i,cost2(:,1),'-o','Color',cm(1,:)); hold on;
for i = 2:ni
    xx = tt(i) + (0:(nIters2 - 1));
    plot(xx,cost2(1:nIters2,i),'-o','Color',cm(2,:));
end
xlabel('time + iteration');
ylabel('cost');
title('onlineUnitaryDil');
axis tight; padAxis();

subplot(2,2,4);
plot(1:nIters2i,repError2(:,1),'-o','Color',cm(1,:)); hold on;
for i = 2:ni
    xx = tt(i) + (0:(nIters2 - 1));
    plot(xx,repError2(1:nIters2,i),'-o','Color',cm(2,:));
end
xlabel('time + iteration');
ylabel('representation error');
title('onlineUnitaryDil');
axis tight; padAxis();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

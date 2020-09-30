%%
% Compare DLS and Online DLS for video denoising
%
% Brian Moore
% brimoor@umich.edu
%

rng(42);

% Data knobs
%path = 'road.mat';
path = 'road396.mat';
idim = [128, nan, 200];             % Raw video dimensions
T    = 5;                           % Streaming size
dt   = 1;                           % Temporal patch stride
pdim = [8, 8, T];                   % Patch sizes
pgap = [2, 2, dt];                  % Patch strides
p    = 0.50;                        % Missing data probability

% Algorithm knobs
lambda   = 0.01;                    % Dictionary regularization parameter
mu       = 0.05;                    % Sparsity regularization
gamma    = 0.9;                     % Forgetting factor (online)
D0       = dctmtx(prod(pdim))';     % Initial dictionary
type     = 'hard';                  % Sparsity regularizer
dr       = 5;                       % Atom rank constraint
ddim     = [prod(pdim(1:2)), T];    % Atom dimensions
unitaryD = false;                   % Unitary dictionary?
fixedD   = false;                   % Fixed dictionary?
flag     = 1;                       % Print stats?
nItersDB = 1;                       % # (D,B) updates per outer iteration
% vvv  dls only  vvv
nIters1  = 10;                      % # outer DLS iters
% vvv  onlineDls only  vvv
nIters2i = nIters1;                 % # outer iters to init online DLS
nIters2  = 10;                      % # outer iters per time for online DLS

% Load video
Xtrue = loadVideo(path,idim,T,dt);
[ny, nx, nt] = size(Xtrue);

% Missing data
M = (rand(ny,nx,nt) > p);
Y = Xtrue;
Y(~M) = 0;

% Interpolated video
Xi = interpVideo(Y,M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts1 = struct();
opts1.A        = 1;
opts1.M        = M;
opts1.xdim     = [ny, nx, nt];
opts1.pdim     = pdim;
opts1.pgap     = pgap;
opts1.type     = type;
opts1.dr       = dr;
opts1.ddim     = ddim;
opts1.unitaryD = unitaryD;
opts1.fixedD   = fixedD;
opts1.nIters   = nIters1;
opts1.nItersDB = nItersDB;
opts1.nItersX  = -1; % Exact updates
opts1.X0       = Xi;
opts1.D0       = D0;
opts1.Xtrue    = Xtrue;
opts1.flag     = flag;
[Xhat1, D1, B1, stats1] = dls(Y,lambda,mu,opts1);
cost1 = stats1.cost(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Online DLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
Xhat2 = zeros(ny,nx,nt);
tt = 1:dt:(nt + 1 - T);
ni = numel(tt);
cost2 = nan(nIters2i,ni);

% Data
M1 = M(:,:,1:T);
Xi1 = Xi(:,:,1:T);
Xtrue1 = Xtrue(:,:,1:T);
Y1 = Y(:,:,1:T);

% Online DLS (t = 1)
opts2 = struct();
opts2.A        = 1;
opts2.M        = M1;
opts2.xdim     = [ny, nx, T];
opts2.pdim     = pdim;
opts2.pgap     = pgap;
opts2.type     = type;
opts2.dr       = dr;
opts2.ddim     = ddim;
opts2.unitaryD = unitaryD;
opts2.fixedD   = fixedD;
opts2.nIters   = nIters2i;
opts2.nItersDB = nItersDB;
opts2.nItersX  = -1; % Exact updates
opts2.X0       = Xi1;
opts2.D0       = D0;
opts2.Xtrue    = Xtrue1;
opts2.accel    = false;
opts2.tau      = 1;
opts2.flag     = flag;
[Xt, Dt, Bt, params, stats2] = onlineDls(Y1,lambda,mu,gamma,opts2);
cost2(:,1) = stats2.cost(:);

% Initialize reconstruction
Xhat2(:,:,1:T) = Xt;

% Online DLS (t = 2,3,...)
for i = 2:ni
    % Data
    t = tt(i);
    Mt = M(:,:,t:(t + T - 1));
    Xit = Xi(:,:,t:(t + T - 1));
    Xtruet = Xtrue(:,:,t:(t + T - 1));
    Yt = Y(:,:,t:(t + T - 1));
    
    % Online DLS
    opts2.M      = Mt;
    opts2.nIters = nIters2;
    %opts2.X0    = zeroFill(Xhat2,Yt,Mt,t,dt);
    opts2.X0     = interpFill(Xhat2,Xit,t,dt);
    opts2.D0     = Dt;
    opts2.B0     = Bt;
    opts2.params = params;
    opts2.Xtrue  = Xtruet;
    [Xt, Dt, Bt, params, stats2] = onlineDls(Yt,lambda,mu,gamma,opts2);
    cost2(1:nIters2,i) = stats2.cost(:);
    
    % Update reconstruction
    Xhat2 = updateRecon(Xhat2,Xt,gamma,t,dt); % gamma-weighted
    %Xhat2 = updateRecon(Xhat2,Xt,0,t,dt);    % latest
    %Xhat2 = updateRecon(Xhat2,Xt,1,t,dt);    % unweighted
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
movie.video = cell2mat({Xtrue,Xi;Xhat1,Xhat2});
opts.xlabels = {'Truth','Interpolated';'Batch','Online'};
PlayMovie(movie,opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cm = linspecer(2);
cfigure();

%{
% DLS
subplot(1,2,1);
plot(1:nIters1,cost1,'-o','Color',cm(1,:)); hold on;
xlabel('iteration');
ylabel('cost');
title('dls');
axis tight; padAxis();

% Online DLS
subplot(1,2,2);
%}
plot(1:nIters2i,cost2(:,1),'-o','Color',cm(1,:)); hold on;
for i = 2:ni
    xx = tt(i) + (0:(nIters2 - 1));
    plot(xx,cost2(1:nIters2,i),'-o','Color',cm(2,:));
end
xlabel('time + iteration');
ylabel('cost');
title('onlineDls');
axis tight; padAxis();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

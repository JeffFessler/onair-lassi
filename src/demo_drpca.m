%%
% Dictionary Robust PCA demo
%
% Brian Moore
% brimoor@umich.edu
%

rng(1);

% Knobs
inpath  = 'mandrill.png';
d       = 40;
lambdaS = 0.1;
lambdaB = 0.01;

% Load data
I      = im2double(imread(inpath));
[m, n] = size(I);

% Generate toy dataset
Y = zeros(m * n,d);
M = false(m * n,d);
for i = 1:d
    % Random observation rectangle
    r = sort(randi(n,1,2));
    c = sort(randi(m,1,2));
    Mi = false(m,n);
    Mi(c(1):c(2),r(1):r(2)) = true;
    
    % Generate frame
    M(:,i) = Mi(:);
    Y(:,i) = I(:) .* M(:,i);
end

% Dictionary Robust PCA
opts.M        = M;          % Frame masks
opts.r        = 1;          % Use OptShrink to generate "robust" panorama
opts.sdim     = [m, n, d];  % Dimensions of S when extracting patches
opts.pdim     = [8, 8, 1];  % Dimensions of a single patch
opts.pgap     = [2, 2, 1];  % Gaps between consecutive patches
opts.type     = 'hard';     % 'hard' or 'soft'
opts.unitaryD = false;      % Unitary dictionary?
opts.nIters   = 1;          % # outer iterations
opts.nItersDB = 1;          % # (D,B) updates per outer iteration
opts.nItersLS = 10;         % # (L,S) updates per outer iteration
opts.D0       = dctmtx(prod(opts.pdim))'; % Initial dictionary
opts.accel    = true;       % Use accelerated proximal gradient?
opts.flag     = 1;          % Print iteration stats?
[L, S, D, B]  = drpca(Y,[],lambdaS,lambdaB,opts);
X = M .* (L + S);

% Visualize result
Y3 = reshape(Y,m,n,d);
X3 = reshape(X,m,n,d);
L3 = reshape(L,m,n,d);
S3 = reshape(S,m,n,d);
opts2.xlabels = {'Y','M .* (L + S)','L','S'};
PlayMovie(cat(2,Y3,X3,L3,S3),opts2);

% Plot dictionary
figure;
dim = [opts.pdim(1), prod(opts.pdim(2:end))];
T0   = tilePatches(opts.D0,dim,dim,1,1,true);
T    = tilePatches(D      ,dim,dim,1,1,true);
subplot(1,2,1); imshow(T0,[]); title('Initial dictionary');
subplot(1,2,2); imshow(T ,[]); title('Final dictionary');

%%
% unitaryDil on real data
%
% Brian Moore
% brimoor@umich.edu
%

rng(42);

% Data knobs
inpath = 'mandrill.png';            % Data
pdim   = [16, 16];                  % Patch sizes
pgap   = [2, 2];                    % Patch strides

% unitaryDil knobs
mu          = 0.5;                  % Sparsity regularization
D0          = dctmtx(prod(pdim))';  % Initial dictionary
opts.type   = 'hard';               % Type of thresholding
opts.D0     = D0;                   % Initial dictionary
opts.nIters = 10;                   % # iterations
opts.flag   = 1;                    % Print iteration status?

% Clustering knobs
kB      = 3;                        % # clusters
nItersB = 50;                       % # k-means iterations
flagB   = true;                     % Print iteration status?

% Load data
I = im2double(imread(inpath));

% Extract patches
[Mp, ~, Np] = patchInds(size(I),pdim,pgap);
Y = I(Mp);

% Learn dictionary
[D, B] = unitaryDil(Y,mu,opts);
T0 = tilePatches(D0,pdim,pdim,1,0,true);
T = tilePatches(D ,pdim,pdim,1,0,true);

% Cluster sparse codes
[~, yB] = kmeans(B,kB,'++',nItersB,flagB);

% Generate cluster assignment map
ZB = clusterMap(size(I),yB,Mp,Np);

% Visualize dictionaries
figure;
subplot(1,2,1); imshow(T0,[]); title('Initial dictionary');
subplot(1,2,2); imshow(T ,[]); title('Final dictionary');

% Visualize cluster memberships
figure;
subplot(1,2,1); imshow(I ,[]); title('Image');
subplot(1,2,2); imshow(ZB,[]); title('Cluster memberships');

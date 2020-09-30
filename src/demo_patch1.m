%% Patch extraction demo

% Knobs
inpath = 'mandrill.png';
pdim = [8, 10];
pgap = pdim; % Non-overlapping

% Load image
X = im2double(imread(inpath));

% Extract patches
adim = size(X);
[Mp, Dp] = patchInds(adim,pdim,pgap);
P = X(Mp);

% Visualize patches
Y = tilePatches(P,pdim,Dp);

% Plot results
figure;
subplot(1,2,1); imshow(X); title('Image');
subplot(1,2,2); imshow(Y); title('Patches');

%% Patch extraction demo (with mask)

% Knobs
inpath = 'mandrill.png';
pdim = [8, 10];
pgap = pdim; % Non-overlapping

% Load image
X = im2double(imread(inpath));

% Mask
[m, n] = size(X);
r = 0.5 * min(m,n);
[XX, YY] = meshgrid(1:n,1:m);
mask = ((XX - 0.5 * n).^2 + (YY - 0.5 * m).^2 < r^2);

% Extract patches
adim = size(X);
[Mp, Dp, Np, Ip] = patchInds(adim,pdim,pgap,mask);
P = zeros(size(Mp,1),prod(Dp));
P(:,Ip) = X(Mp);

% Visualize patches
Y = tilePatches(P,pdim,Dp);

% Plot results
figure;
subplot(1,3,1); imshow(X); title('Image');
subplot(1,3,2); imshow(mask); title('Mask');
subplot(1,3,3); imshow(Y); title('Patches');

%% Speed test (vs. im2col)

% Knobs
inpath = 'mandrill.png';
pdim   = [8, 10];

% Load image
X = im2double(imread(inpath));

%--------------------------------------------------------------------------
% Distinct test
%--------------------------------------------------------------------------
disp('Distinct test');

% Method #1
tic;
M1 = patchInds(size(X),pdim,pdim);
P1 = X(M1);
patchIndsTime = toc %#ok

% Method #2
tic;
P2 = im2col(X,pdim,'distinct');
im2colTime = toc %#ok

% Compute difference
error = norm(P1 - P2) %#ok
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sliding test
%--------------------------------------------------------------------------
disp('Sliding test');

% Method #1
tic;
M1 = patchInds(size(X),pdim,[1 1]);
P1 = X(M1);
patchIndsTime = toc %#ok

% Method #2
tic;
P2 = im2col(X,pdim,'sliding');
im2colTime = toc %#ok

% Compute difference
error = norm(P1 - P2) %#ok
%--------------------------------------------------------------------------

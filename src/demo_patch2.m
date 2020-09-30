%% Patch-based compression demo
%

%% Spot check

rng(42);

% Knobs
inpath = 'cameraman.png';   % Input image
pdim = [8, 8];              % Patch dimensions
pgap = [1, 1];              % Patch strides
r1 = 15;                    % All-patch rank, <= prod(pdim)                              
r2 = 1;                     % Per-patch rank, <= min(pdim)
SNR = 15;                   % SNR, in dB

% Load image
Xtrue = im2double(imread(inpath));
xdim = size(Xtrue);

% Add noise
Y = corrupt(Xtrue, SNR);

% Extract patches
[Mp, Dp, Np] = patchInds(xdim, pdim, pgap);
P = Y(Mp);

% All-patch TSVD
P1 = tsvd(P, r1);

% Per-patch TSVDs
P2 = zeros(size(P));
for i = 1:size(P2,2)
    pr2i = tsvd(reshape(P(:,i), pdim), r2);
    P2(:,i) = pr2i(:);
end

% Reconstruct images
Xhat1 = patchRecon(P1, Mp, Np);
Xhat2 = patchRecon(P2, Mp, Np);

% Compute errors
pSNRfcn = @(Xhat, Xtrue) 10 * log10(max(abs(Xtrue(:)))^2 / ...
                                    mean(abs(Xhat(:) - Xtrue(:)).^2));
pSNRY = pSNRfcn(Y, Xtrue);
pSNR1 = pSNRfcn(Xhat1, Xtrue);
pSNR2 = pSNRfcn(Xhat2, Xtrue);

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
%{
% Visualize patches
figure;
imshow(tilePatches(P, pdim, Dp, [], [] , false), [0, 1]);
title('Patches');
%}

figure;

% True image
subplot(2,2,1);
imshow(Xtrue, [0, 1]);
title('True image');

% Noisy image
subplot(2,2,2);
imshow(Y, [0, 1]);
title(sprintf('Noisy image [pSNR = %.1fdB]', pSNRY));

% All-patch TSVD
subplot(2,2,3);
imshow(Xhat1, [0, 1]);
title(sprintf('All-patch TSVD [pSNR = %.1fdB]', pSNR1));

% Per-patch TSVDs
subplot(2,2,4);
imshow(Xhat2, [0, 1]);
title(sprintf('Per-patch TSVDs [pSNR = %.1fdB]', pSNR2));

%{
% Save figure
export_fig('-pdf', '-transparent', 'demo2_recons');
%}
%--------------------------------------------------------------------------

%% Parameter sweep

% Knobs
inpath = 'cameraman.png';               % Input image
pdim = [8, 8];                          % Patch dimensions
pgap = [1, 1];                          % Patch strides
r1 = [1, 3, 6, 9, 12, 15, 20, 40, 64];  % All-patch rank, <= prod(pdim)                              
r2 = 1:8;                               % Per-patch rank, <= min(pdim)
SNR = [5, 10, 15];                      % SNR, in dB

% Load image
Xtrue = im2double(imread(inpath));
xdim = size(Xtrue);

% Generate noise
N = randn(xdim);

% PNR function
pSNRfcn = @(Xhat, Xtrue) 10 * log10(max(abs(Xtrue(:)))^2 / ...
                                    mean(abs(Xhat(:) - Xtrue(:)).^2));

% Perform experiments
nr1 = numel(r1);
nr2 = numel(r2);
nSNR = numel(SNR);
pSNR1 = nan(nr1,nSNR);
pSNR2 = nan(nr2,nSNR);
for j = 1:nSNR
    % Add noise
    sigma = 10^(-SNR(j) / 20) * norm(Xtrue(:)) / sqrt(numel(Xtrue));
    Y = Xtrue + sigma * N;
    
    % Extract patches
    [Mp, Dp, Np] = patchInds(xdim, pdim, pgap);
    P = Y(Mp);
    
    % All-patch TSVDs
    for i = 1:nr1
        P1 = tsvd(P, r1(i));
        Xhat1 = patchRecon(P1, Mp, Np);
        pSNR1(i, j) = pSNRfcn(Xhat1, Xtrue);
    end
    
    % Per-patch TSVDs
    for i = 1:nr2
        P2 = zeros(size(P));
        for k = 1:size(P2,2)
            pr2k = tsvd(reshape(P(:,k), pdim), r2(i));
            P2(:,k) = pr2k(:);
        end
        Xhat2 = patchRecon(P2, Mp, Np);
        pSNR2(i, j) = pSNRfcn(Xhat2, Xtrue);
    end
end

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
%cfigure([279, 227]);
cfigure([900, 235]);
cm = linspecer(2);
ploth = @semilogx;

% Plot results
phndl = zeros(1, 2);
[nh, nw] = bestSubplotShape(nSNR);
for j = 1:nSNR
    subplot(nh, nw, j);
    phndl(1) = ploth(r1, pSNR1(:,j), '-o', 'Color', cm(1, :)); hold on;
    phndl(2) = ploth(r2, pSNR2(:,j), '-x', 'Color', cm(2, :));
    title(sprintf('SNR = %.0fdB', SNR(j)));
    xlabel('rank (r)');
    ylabel('pSNR, in dB');
    axis tight; padAxis();
end
legend(phndl, 'All-patch TSVD', 'Per-patch TSVDs', 'Location', 'SE');
matchAxes([], 0, 'y');

%{
% Save figure
export_fig('-pdf', '-transparent', 'demo2_psnrs');
%}
%--------------------------------------------------------------------------

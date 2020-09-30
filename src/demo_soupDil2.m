%%
% soupDil on synthetic data
%
% Brian Moore
% brimoor@umich.edu
%

rng(42);

% Data knobs
d   = 64;
m   = 64;
n   = 4096;
p   = 0.10;
SNR = 0;

% soupDil knobs
mu     = 0.6;           % Sparsity regularization
nIters = 50;            % # iterations
type   = 'hard';        % Type of thresholding
flag   = 0;             % Print iteration status?

% Generate data
Dtrue = dctmtx(d)';
Btrue = randn(m,n) .* (rand(m,n) < p);
Y = corrupt(Dtrue * Btrue,SNR);

% Initial dictionary
D0 = corrupt(Dtrue,SNR);
D0 = bsxfun(@rdivide,D0,sqrt(sum(abs(D0).^2,1))); % unit-norm columns

% SOUP-DIL
opts = struct();
opts.nIters = nIters;
opts.D0     = D0;
opts.type   = type;
opts.tau    = 0.0175; % Sparsity optimization
opts.flag   = flag;
[D, B, stats] = soupDil(Y,mu,opts);

% Errors
NRMSEfcn = @(X,Xtrue) norm(X(:) - Xtrue(:)) / norm(Xtrue(:));
nrmseD = NRMSEfcn(D,Dtrue) %#ok
nrmseB = NRMSEfcn(B,Btrue) %#ok

%--------------------------------------------------------------------------
% Plot results
%--------------------------------------------------------------------------
cfigure();
cm = linspecer(2);

% cost
subplot(2,2,1);
plot(1:nIters,stats.cost,'b-o');
xlabel('iteration');
title('cost');
axis tight; padAxis();

% sparsity
subplot(2,2,2);
plot(1:nIters,stats.sparsity,'b-o');
xlabel('iteration');
title('sparsity');
axis tight; padAxis();

% deltaD
subplot(2,2,3);
semilogy(1:nIters,stats.deltaD,'b-o');
xlabel('iteration');
title('deltaD');
axis tight; padAxis();

% deltaB
subplot(2,2,4);
semilogy(1:nIters,stats.deltaB,'b-o');
xlabel('iteration');
title('deltaB');
axis tight; padAxis();
%--------------------------------------------------------------------------

function [L, S, sL, stats] = patchRPCA(Y,Z,lambdaL,lambdaS,Mp,Np,varargin)
%
% Syntax:       [L, S] = patchRPCA(Y,Z,lambdaL,lambdaS,Mp,Np);
%               [L, S] = patchRPCA(Y,Z,lambdaL,lambdaS,Mp,Np,opts);
%               [L, S, sL, stats] = patchRPCA(Y,Z,lambdaL,lambdaS,Mp,Np);
%               [L, S, sL, stats] = patchRPCA(Y,Z,lambdaL,lambdaS,Mp,Np,opts);
%               
% Inputs:       Y is an m x n data matrix
%               
%               Z is a b x d data matrix
%               
%               lambdaL >= 0 is the nuclear norm regularization parameter.
%               Note that lambdaL is ignored when ~isnan(r)
%               
%               lambdaS >= 0 is the patch-sparsity regularization parameter
%               
%               Mp is a b x d patch index matrix from patchInds()
%               
%               Np is a p x n patch count matrix from patchInds()
%               
%               [OPTIONAL] opts is a struct containing one or more of the
%               following fields. The default values are in ()
%                   
%                   opts.A (1) is an m x p system matrix
%                   
%                   opts.M (1) is an m x n weight matrix
%                   
%                   opts.r (nan) is the OptShrink rank parameter. When a
%                   non-nan r is specified, lambdaL is ignored and
%                   OptShrink is used in place of SVT for all L updates
%                   
%                   opts.nIters (50) is the number of iterations to perform
%                   
%                   opts.L0 (A' * Y) is a p x n matrix specifying the
%                   initial L iterate
%                   
%                   opts.S0 (zeros(p,n)) is a p x n matrix specifying the
%                   initial S iterate
%                   
%                   opts.Xtrue (nan) is the ground truth X = L + S matrix
%                   to use for MSE calculations
%                   
%                   opts.NRMSEfcn (@computeNRMSE) is the function to use
%                   when computing the NRMSE of X = L + S after each
%                   iteration
%                   
%                   opts.accel (true) specifies whether to use accelerated
%                   proximal gradient steps
%                   
%                   opts.tau 0.99 / ((1 + accel) * norm(A)^2) is the step
%                   size parameter, and should satisfy
%                   
%                       tau <= 1 / (2 * norm(A)^2), when accel = true
%                       tau <  1 /      norm(A)^2 , when accel = false
%                   
%                   opts.flag (1) determines what status updates to print
%                   to the command window. The choices are
%                   
%                       flag = 0: no printing
%                       flag = 1: print iteration status
%               
% Outputs:      L is a p x n matrix containing the estimated low-rank
%               component
%               
%               S is a p x n matrix containing the estimated patch-sparse
%               component
%               
%               sL is a min(p,n) x 1 vector containing the singular values
%               of L
%               
%               stats is a struct containing the following fields:
%               
%                   stats.nIters is the number of iterations performed
%                   
%                   stats.cost is a 1 x nIters vector containing the cost
%                   function values at each iteration
%                   
%                   stats.nrmse is a 1 x nIters vector containing the NRMSE
%                   of X = L + S with respect to Xtrue at each iteration
%                   
%                   stats.delta is a 1 x nIters vector containing the
%                   relative convergence of X = L + S at each iteration,
%                   defined as \|X_{k + 1} - X_k\|_F / \|X_k\|_F
%                   
%                   stats.time is a 1 x nIters vector containing the time,
%                   in seconds, required to perform each outer iteration
%               
% Description:  Solves the Patch Robust PCA (PRPCA) problem
%               
%               min_{L,S} 0.5\|\sqrt{M} \odot (Y - A(L + S))\|_F^2 + 
%                            \lambda_L\|L\|_{\star} + 
%                         0.5\lambda_S\|P(S) - Z\|_F^2
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         April 8, 2016
%               June 6, 2017
%

% Parse inputs
[A, M, r, nIters, L, S, Xtrue, NRMSEfcn, accel, tau, flag] = ...
                                                parseInputs(Y,varargin{:});
PRINT_STATS   = (flag > 0);
COMPUTE_STATS = PRINT_STATS || (nargout == 4);
USE_OPTSHRINK = ~isnan(r);

% Initialize stats
if COMPUTE_STATS
    % Cost function
    Psi = @(X,sL,S) 0.5 * norm(vec(sqrt(M) .* (Y - A * X)))^2 + ...
                lambdaL * sum(sL) + 0.5 * lambdaS * norm(vec(S(Mp) - Z))^2;
    
    % Stats-printing function
    iterFmt = sprintf('%%0%dd',ceil(log10(nIters + 1)));
    out     = printFcn('Iter' ,iterFmt,'cost','%.2f','nrmse','%.3f', ...
                       'delta','%.3e' ,'time','%.2fs');
    
    % Initialize stats
    X     = L + S;
    cost  = nan(1,nIters);
    nrmse = nan(1,nIters);
    delta = nan(1,nIters);
    time  = nan(1,nIters);
end

% Construct patch approximation
PiZ = patchApprox(Mp,Z,size(S));

% Patch Robust PCA
if PRINT_STATS
    fprintf('***** Patch Robust PCA *****\n');
end
if accel
    % Initialize accelerated method
    t     = 0;
    Llast = L;
    Slast = S;
end
for it = 1:nIters
    % Initilize iteration
    itimer = tic;
    
    % Proximal gradient update
    if accel
        % Accelerated step
        tlast = t;
        t     = 0.5 * (1 + sqrt(1 + 4 * t^2));
        Lbar  = L + ((tlast - 1) / t) * (L - Llast);
        Sbar  = S + ((tlast - 1) / t) * (S - Slast);
        Llast = L;
        Slast = S;
        Gbar  = A' * (M .* (A * (Lbar + Sbar) - Y));
        if USE_OPTSHRINK
            % OptShrink
            [L, sL] = OptShrink(Lbar - tau * Gbar,r);
        else
            % SVT
            [L, sL] = SVT(Lbar - tau * Gbar,tau * lambdaL);
        end
        S = patchProx(Sbar - tau * Gbar,tau * lambdaS,PiZ,Np);
    else
        % Standard step
        G = A' * (M .* (A * (L + S) - Y));
        if USE_OPTSHRINK
            % OptShrink
            [L, sL] = OptShrink(L - tau * G,r);
        else
            % SVT
            [L, sL] = SVT(L - tau * G,tau * lambdaL);
        end
        S = patchProx(S - tau * G,tau * lambdaS,PiZ,Np);
    end
    
    % Record stats
    if COMPUTE_STATS
        Xlast = X;
        X     = L + S;
        if ~USE_OPTSHRINK
            cost(it) = Psi(X,sL,S);
        end
        nrmse(it) = NRMSEfcn(X,Xtrue);
        delta(it) = computeNRMSE(X,Xlast);
        time(it)  = toc(itimer);
        if PRINT_STATS
            out(it,cost(it),nrmse(it),delta(it),time(it)); 
        end
    end
end

% Return stats
if COMPUTE_STATS
    stats.nIters = nIters;
    stats.cost   = cost;
    stats.nrmse  = nrmse;
    stats.delta  = delta;
    stats.time   = time;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Singular value thresholding
function [Y, sY] = SVT(X,lambda)
[UX, SX, VX] = svd(X,'econ');
sY = soft(diag(SX),lambda);
Y  = UX * diag(sY) * VX';

% Soft thresholding
function Y = soft(X,lambda)
Y = sign(X) .* max(abs(X) - lambda,0);

% Patch proximal operator
% \min_{X} \|Y - X\|_F^2 + 0.5 \lambda \|P(X) - Z\|_F^2
function X = patchProx(Y,lambda,PiZ,Np)
X = (Y + lambda * PiZ) ./ (1 + lambda * Np);

% Patch approximation matrix
% PiZ = \sum_{j=1}^n P_j^T z_j
function PiZ = patchApprox(Mp,Z,dim)
PiZ = zeros(dim);
n = size(Mp,2);
for i = 1:n
    PiZ(Mp(:,i)) = PiZ(Mp(:,i)) + Z(:,i);
end

% Print function
function out = printFcn(varargin)
str = [sprintf('%s[%s] ',varargin{:}), '\n'];
out = @(varargin) fprintf(str,varargin{:});

% Compute NRMSE
function err = computeNRMSE(Xhat,X)
denom = norm(X(:));
if isnan(denom)
    err = nan;
elseif denom == 0
    err = 0;
else
    err = norm(Xhat(:) - X(:)) / denom;
end

% Vectorize input
function x = vec(X)
x = X(:);

% Parse inputs
function [A, M, r, nIters, L0, S0, Xtrue, NRMSEfcn, accel, tau, flag] = ...
                                                        parseInputs(Y,opts)
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

% Standard args
A        = parseField(opts,'A',1);
M        = parseField(opts,'M',1);
r        = parseField(opts,'r',nan);
nIters   = parseField(opts,'nIters',50);
L0       = parseField(opts,'L0',nan);
S0       = parseField(opts,'S0',nan);
Xtrue    = parseField(opts,'Xtrue',nan);
NRMSEfcn = parseField(opts,'NRMSEfcn',@computeNRMSE);
accel    = parseField(opts,'accel',true);
tau      = parseField(opts,'tau',nan);
flag     = parseField(opts,'flag',1);

% Expensive defaults
if isnan(L0),  L0   = A' * Y;                           end
if isnan(S0),  S0   = zeros(size(L0));                  end
if isnan(tau), tau  = 1 / ((1 + accel) * norm(A)^2);    end

% Parse struct field
function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end

function [X, stats] = patchLS(Y,Z,lambda,Mp,Np,varargin)
%
% Syntax:       X = patchLS(Y,Z,lambda,Mp,Np);
%               X = patchLS(Y,Z,lambda,Mp,Np,opts);
%               [X, stats] = patchLS(Y,Z,lambda,Mp,Np);
%               [X, stats] = patchLS(Y,Z,lambda,Mp,Np,opts);
%               
% Inputs:       Y is an m x n data matrix
%               
%               Z is a b x d data matrix
%               
%               lambda >= 0 is the patch regularization parameter
%               
%               Mp is an b x d patch index matrix from patchInds()
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
%                   opts.nIters (50) is the number of iterations to
%                   perform. Note that, when A = 1, the exact solution is
%                   always computed, so nIters is ignored
%                   
%                   opts.X0 (A' * Y) is a p x n matrix specifying the
%                   initial X iterate
%                   
%                   opts.Xtrue (nan) is the ground truth X matrix to use 
%                   for MSE calculations
%                   
%                   opts.NRMSEfcn (@computeNRMSE) is the function to use
%                   when computing the NRMSE of X after each iteration
%                   
%                   opts.accel (true) specifies whether to use accelerated
%                   proximal gradient steps
%                   
%                   opts.tau ((0.99 + ~accel) / norm(A)^2) is the step size
%                   parameter, and should satisfy
%                   
%                       tau <= 1 / norm(A)^2, when accel = true
%                       tau <  2 / norm(A)^2, when accel = false
%                   
%                   opts.flag (1) determines what status updates to print
%                   to the command window. The choices are
%                   
%                       flag = 0: no printing
%                       flag = 1: print iteration status
%               
% Outputs:      X is the estimated p x n matrix
%               
%               stats is a struct containing the following fields:
%               
%                   stats.nIters is the number of iterations performed
%                   
%                   stats.cost is a 1 x nIters vector containing the cost
%                   function values at each iteration
%                   
%                   stats.nrmse is a 1 x nIters vector containing the NRMSE
%                   of X with respect to Xtrue at each iteration
%                   
%                   stats.delta is a 1 x nIters vector containing the
%                   relative convergence of X at each iteration, defined as
%                   \|X_{k + 1} - X_k\|_F / \|X_k\|_F
%                   
%                   stats.time is a 1 x nIters vector containing the time,
%                   in seconds, required to perform each outer iteration
%               
% Description:  Solves the Patch Least-Squares (PLS) problem
%               
%               \min_{X} 0.5 \|\sqrt{M} \odot (Y - AX)\|_F^2 +
%                        0.5 \lambda \|P(X) - Z\|_F^2
%               
% Note:         When A = 1 (inpainting), the exact solution is computed
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         April 8, 2016
%               January 24, 2017
%               January 31, 2017
%               June 6, 2017
%

% Parse inputs
[A, M, nIters, X, Xtrue, NRMSEfcn, accel, tau, flag] = ...
                                                parseInputs(Y,varargin{:});
PRINT_STATS   = (flag > 0);
COMPUTE_STATS = PRINT_STATS || (nargout == 2);

% Initialize stats
if COMPUTE_STATS
    % Cost function
    Psi = @(X) 0.5 * norm(vec(sqrt(M) .* (Y - A * X)))^2 + ...
               0.5 * lambda * norm(vec(X(Mp) - Z))^2;
    
    % Stats-printing function
    iterFmt = sprintf('%%0%dd',ceil(log10(nIters + 1)));
    out     = printFcn('Iter' ,iterFmt,'cost','%.2f','nrmse','%.3f', ...
                       'delta','%.3e' ,'time','%.2fs');
    
    % Initialize stats
    cost  = nan(1,nIters);
    nrmse = nan(1,nIters);
    delta = nan(1,nIters);
    time  = nan(1,nIters);
end

% Construct patch approximation
PiZ = patchApprox(Mp,Z,size(X));

% Patch Least-Squares
if PRINT_STATS
    fprintf('***** Patch Least-Squares *****\n');
end
if isequal(A,1)
    % Handle inpainting (A == 1) separately
    if nIters > 0
        % Warn user that nIters > 0 is unnecessary
        warning('nIters > 0 is unnecessary when A == 1'); %#ok
    end
    
    % Exact solution
    itimer = tic;
    X = weightedPatchProx(Y,M,lambda,PiZ,Np);
    
    % Record stats
    if COMPUTE_STATS
        cost  = Psi(X);
        nrmse = NRMSEfcn(X,Xtrue);
        delta = 0;
        time  = toc(itimer);
        if PRINT_STATS
            out(0,cost,nrmse,delta,time); 
        end
    end
else
    % Handle general case (A ~= 1) via proximal gradient
    if accel
        % Initialize accelerated method
        t     = 0;
        Xlast = X;
    end
    for it = 1:nIters
        % Initilize iteration
        itimer = tic;
        
        % Proximal gradient update
        if accel
            % Accelerated step
            tlast = t;
            t     = 0.5 * (1 + sqrt(1 + 4 * t^2));
            Xbar  = X + ((tlast - 1) / t) * (X - Xlast);
            Xlast = X;
            Gbar  = A' * (M .* (A * Xbar - Y));
            X     = patchProx(Xbar - tau * Gbar,tau * lambda,PiZ,Np);
        else
            % Standard step
            Xlast = X;
            G     = A' * (M .* (A * X - Y));
            X     = patchProx(X - tau * G,tau * lambda,PiZ,Np);
        end
        
        % Record stats
        if COMPUTE_STATS
            cost(it)  = Psi(X);
            nrmse(it) = NRMSEfcn(X,Xtrue);
            delta(it) = computeNRMSE(X,Xlast);
            time(it)  = toc(itimer);
            if PRINT_STATS
                out(it,cost(it),nrmse(it),delta(it),time(it)); 
            end
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

% Patch proximal operator
% \min_{X} \|Y - X\|_F^2 + 0.5 \lambda \|P(X) - Z\|_F^2
function X = patchProx(Y,lambda,PiZ,Np)
X = (Y + lambda * PiZ) ./ (1 + lambda * Np);

% Weighted patch proximal operator
% \min_{X} \|\sqrt{M} \odot (Y - X)\|_F^2 + \lambda\|P(X) - Z\|_F^2
function X = weightedPatchProx(Y,M,lambda,PiZ,Np)
X = (M .* Y + lambda * PiZ) ./ (M + lambda * Np);

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
function [A, M, nIters, X0, Xtrue, NRMSEfcn, accel, tau, flag] = ...
                                                        parseInputs(Y,opts)
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

% Standard args
A        = parseField(opts,'A',1);
M        = parseField(opts,'M',1);
nIters   = parseField(opts,'nIters',50);
X0       = parseField(opts,'X0',nan);
Xtrue    = parseField(opts,'Xtrue',nan);
NRMSEfcn = parseField(opts,'NRMSEfcn',@computeNRMSE);
accel    = parseField(opts,'accel',true);
tau      = parseField(opts,'tau',nan);
flag     = parseField(opts,'flag',1);

% Expensive defaults
if isnan(X0),  X0   = A' * Y;                      end
if isnan(tau), tau  = (0.99 + ~accel) / norm(A)^2; end

% Parse struct field
function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end

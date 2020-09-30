function [D, B, stats] = unitaryDil(Y,mu,varargin)
%
% Syntax:       [D, B] = unitaryDil(Y,mu);
%               [D, B] = unitaryDil(Y,mu,opts);
%               [D, B, stats] = unitaryDil(Y,mu);
%               [D, B, stats] = unitaryDil(Y,mu,opts);
%               
% Inputs:       Y is a d x n data matrix
%               
%               mu >= 0 is the sparsity regularization parameter
%               
%               [OPTIONAL] opts is a struct containing one or more of the
%               following fields. The default values are in ()
%                   
%                   opts.D0 (dctmtx(d)') is a d x d matrix containing the
%                   initial unitary dictionary
%                   
%                   opts.B0 (zeros(d,n)) is a d x n matrix containing the
%                   initial sparse codes
%                   
%                   opts.type ('hard') can be 'hard' or 'soft' and
%                   specifies whether apply hard or soft thresholding to B
%                   
%                   opts.fixedD (false) determines whether to fix the
%                   initial dictionary (true) or learn it (false)
%                   
%                   opts.nIters (50) is the number of block coordinate
%                   descent iterations to perform
%                   
%                   opts.flag (1) determines what status updates to print
%                   to the command window. The choices are
%                   
%                       flag = 0: no printing
%                       flag = 1: print iteration status
%               
% Outputs:      D is a d x d unitary matrix containing the final dictionary
%               
%               B is an d x n matrix containing the final sparse codes
%               
%               stats is a struct containing the following fields:
%               
%                   stats.nIters is the number of iterations performed
%                   
%                   stats.cost is a 1 x nIters vector containing the cost
%                   function values at each iteration
%                   
%                   stats.repError is a 1 x nIters vector containing the
%                   representation error at each iteration, defined as
%                   \|Y - D_k B_k\|_F / \|Y\|_F
%                   
%                   stats.deltaD is a 1 x nIters vector containing the
%                   relative convergence of D at each iteration, defined as
%                   \|D_k - D_{k - 1}\|_F / \|D_{k - 1}\|_F
%                   
%                   stats.deltaB is a 1 x nIters vector containing the
%                   relative convergence of B at each iteration, defined as
%                   \|B_k - B_{k - 1}\|_F / \|B_{k - 1}\|_F
%                   
%                   stats.sparsity is a 1 x nIters vector containing the
%                   sparsity, in percent, of B at each iteration
%                   
%                   stats.time is a 1 x nIters vector containing the time,
%                   in seconds, required to perform each outer iteration
%               
% Description:  Performs dictionary learning via the Unitary DIctionary
%               Learning (Unitary-DIL) problem. That is,
%               when type == 'hard':
%               
%               min_{D,B}  0.5 \|Y - D B\|_F^2 + 0.5 \mu^2 \|B\|_0
%               
%               subject to D' * D == I
%               
%               Or, when type == 'soft':
%               
%               min_{D,B}  0.5 \|Y - D B\|_F^2 + \mu \|B\|_1
%               
%               subject to D' * D == I
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         May 15, 2017
%

% Parse inputs
[D, B, type, fixedD, nIters, flag] = parseInputs(Y,varargin{:});
PRINT_STATS     = (flag > 0);
COMPUTE_STATS   = PRINT_STATS || (nargout == 3);
USE_HARD_THRESH = strcmpi(type,'hard');
OPTIMIZE_DICT   = ~fixedD;

% Initialize stats
if COMPUTE_STATS
    % Cost function
    Yfit = @(D,B) 0.5 * norm(vec(Y - D * B))^2;
    if strcmpi(type,'hard')
        % Ell-0 regularization
        Bfit = @(B) 0.5 * mu^2 * nnz(B);
    elseif strcmpi(type,'soft')
        % Ell-1 regularization
        Bfit = @(B) mu * norm(vec(B),1);
    else
        % Unsupported shrinkage
        error('Unsupported shrinkage type ''%s''',type);
    end
    Psi = @(D,B) Yfit(D,B) + Bfit(B);
    
    % Stats-printing function
    iterFmt = sprintf('%%0%dd',ceil(log10(nIters + 1)));
    out     = printFcn('Iter'    ,iterFmt,'cost'    ,'%.2f'  , ...
                       'repError','%.3f' ,'deltaD'  ,'%.3e'  , ...
                       'deltaB'  ,'%.3e' ,'sparsity','%.2f%%', ...
                       'time'    ,'%.2fs');
    
    % Initialize stats
    cost     = nan(1,nIters);
    repError = nan(1,nIters);
    deltaD   = nan(1,nIters);
    deltaB   = nan(1,nIters);
    sparsity = nan(1,nIters);
    time     = nan(1,nIters);
end

% Unitary-DIL
if PRINT_STATS
    fprintf('***** Unitary-DIL *****\n');
end
for it = 1:nIters
    % Initialize iteration
    itimer = tic;
    if COMPUTE_STATS
        % Save last iterates
        Blast = B;
        Dlast = D;
    end
    
    % B-update
    if USE_HARD_THRESH
        % Hard thresholding
        B = hard(D' * Y,mu);
    else
        % Soft thresholding
        B = soft(D' * Y,mu);
    end
    
    % D-update
    if OPTIMIZE_DICT
        D = procrustes(Y,B);
    end
    
    % Record stats
    if COMPUTE_STATS
        cost(it)     = Psi(D,B);
        repError(it) = computeNRMSE(D * B,Y);
        deltaD(it)   = computeNRMSE(D,Dlast);
        deltaB(it)   = computeNRMSE(B,Blast);
        sparsity(it) = 100 * (nnz(B) / numel(B));
        time(it)     = toc(itimer);
        if PRINT_STATS
            out(it,cost(it),repError(it),deltaD(it),deltaB(it), ...
                sparsity(it),time(it));
        end
    end
end

% Return stats
if COMPUTE_STATS
    stats.nIters   = nIters;
    stats.cost     = cost;
    stats.repError = repError;
    stats.deltaD   = deltaD;
    stats.deltaB   = deltaB;
    stats.sparsity = sparsity;
    stats.time     = time;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hard thresholding
function Y = hard(X,mu)
Y = X .* (abs(X) > mu);

% Soft thresholding
function Y = soft(X,mu)
Y = sign(X) .* max(abs(X) - mu,0);

% Procrustes solution to
% min_{Q: Q'Q = I} \|X - Q Y\|_F
function Q = procrustes(X,Y)
[U, ~, V] = svd(X * Y','econ');
Q = U * V';

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
function [D0, B0, type, fixedD, nIters, flag] = parseInputs(Y,opts)
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

% Standard args
D0     = parseField(opts,'D0',nan);
B0     = parseField(opts,'B0',nan);
type   = parseField(opts,'type','hard');
fixedD = parseField(opts,'fixedD',false);
nIters = parseField(opts,'nIters',50);
flag   = parseField(opts,'flag',1);

% Expensive defaults
if isnan(D0), D0 = dctmtx(size(Y,1))';          end
if isnan(B0), B0 = zeros(size(D0,2),size(Y,2)); end

% Parse struct field
function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end

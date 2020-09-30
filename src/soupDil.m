function [D, B, stats] = soupDil(Y,mu,varargin)
%
% Syntax:       [D, B] = soupDil(Y,mu);
%               [D, B] = soupDil(Y,mu,opts);
%               [D, B, stats] = soupDil(Y,mu);
%               [D, B, stats] = soupDil(Y,mu,opts);
%               
% Inputs:       Y is a d x n data matrix
%               
%               mu >= 0 is the sparsity regularization parameter
%               
%               [OPTIONAL] opts is a struct containing one or more of the
%               following fields. The default values are in ()
%                   
%                   opts.D0 (dctmtx(d)') is a d x m matrix containing the
%                   initial dictionary. Overcomplete dictionaries (m > d) 
%                   are allowed
%                   
%                   opts.B0 (zeros(m,n)) is an m x n matrix containing the
%                   initial sparse codes
%                   
%                   opts.type ('hard') can be 'hard' or 'soft' and
%                   specifies whether apply hard or soft thresholding to B
%                   
%                   opts.dr (nan) is a rank constraint on the dictionary
%                   atoms. Note: both dr and ddim must be non-nan to apply
%                   a rank constraint
%                   
%                   opts.ddim (nan) is a 1 x 2 vector describing how to
%                   reshape the dictionary atoms into a matrix before 
%                   applying the rank constraint. Note: both dr and ddim
%                   must be non-nan to apply a rank constraint
%                   
%                   opts.fixedD (false) determines whether to fix the
%                   initial dictionary (true) or learn it (false)
%                   
%                   opts.nIters (50) is the number of block coordinate
%                   descent iterations to perform
%                   
%                   opts.order (1:m) is a 1 x m vector specifying the order
%                   in which the columns of D and rows of B are updated. 
%                   Set order = -1 to randomize the order each iteration
%                   
%                   opts.tau (-1) is the sparsity threshold to use to
%                   optimize the storage of B. By default, B is a full
%                   matrix, but, when tau > 0, B is converted to a sparse
%                   matrix for iterations where nnz(B) / numel(B) > tau.
%                   Empirically, a good value to try is tau = 0.015
%                   
%                   opts.flag (1) determines what status updates to print
%                   to the command window. The choices are
%                   
%                       flag = 0: no printing
%                       flag = 1: print iteration status
%               
% Outputs:      D is a d x m matrix containing the final dictionary
%               
%               B is an m x n matrix containing the final sparse codes
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
% Description:  Performs dictionary learning via the Sum of OUter Products
%               DIctionary Learning (SOUP-DIL) problem. That is, when
%               type == 'hard':
%               
%               min_{D,B}  0.5 \|Y - D B\|_F^2 + 0.5 \mu^2 \|B\|_0
%               
%               subject to \|D(:,k)\|_2 = 1
%                          rank(R(D(:,k)) <= dr
%               
%               Or, when type == 'soft':
%               
%               min_{D,B}  0.5 \|Y - D B\|_F^2 + \mu \|B\|_1
%               
%               subject to \|D(:,k)\|_2 = 1
%                          rank(R(D(:,k)) <= dr
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         June 17, 2016
%               November 6, 2016
%               January 25, 2017
%               January 31, 2017
%               May 15, 2017
%

% Parse inputs
[D, B, type, dr, ddim, fixedD, nIters, order, tau, flag] = ...
                                                parseInputs(Y,varargin{:});
RANDOMIZE_ORDER  = isequal(order,-1);
PRINT_STATS      = (flag > 0);
COMPUTE_STATS    = PRINT_STATS || (nargout == 3);
USE_HARD_THRESH  = strcmpi(type,'hard');
OPTIMIZE_STORAGE = (tau > 0);
LOW_RANK_ATOMS   = (dr < min(ddim));
OPTIMIZE_DICT    = ~fixedD;

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

% SOUP-DIL
if PRINT_STATS
    fprintf('***** SOUP-DIL *****\n');
end
[d, m] = size(D);
d0 = [1; zeros(d - 1,1)]; % Zero-reset atom
for it = 1:nIters
    % Initialize iteration
    itimer = tic;
    if RANDOMIZE_ORDER
        % Randomize order
        order = randperm(m);
    end
    if OPTIMIZE_STORAGE
        % Optimize for sparsity
        B = optimizeStorage(B,tau);
    end
    if COMPUTE_STATS
        % Save last iterates
        Blast = B;
        Dlast = D;
    end
    
    % Block coordinate descent
    DtY = D' * Y; % Precompute for efficiency
    for k = order        
        % Update B(k,:)
        bk = DtY(k,:) - (D(:,k)' * D) * B + B(k,:);
        if USE_HARD_THRESH
            % Hard thresholding
            Mk = (abs(bk) > mu);
            bk(~Mk) = 0;
        else
            % Soft thresholding
            bk = sign(bk) .* max(abs(bk) - mu,0);
            Mk = (bk ~= 0);            
        end
        B(k,:) = bk;
        ANY_MK = any(Mk);
        
        % Update D(:,k)
        if OPTIMIZE_DICT
            if ANY_MK
                bMk = bk(Mk)';
                BbMk = B(:,Mk) * bMk;
                dk = Y(:,Mk) * bMk - D * BbMk + BbMk(k) * D(:,k);
            else
                dk = d0;
            end
            if LOW_RANK_ATOMS
                % Low-rank atoms
                [Uk, Sk, Vk] = svd(reshape(dk,ddim),'econ');
                dk = vec(Uk(:,1:dr) * Sk(1:dr,1:dr) * Vk(:,1:dr)');
            end
            D(:,k) = dk / norm(dk);
        end
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

% Print function
function out = printFcn(varargin)
str = [sprintf('%s[%s] ',varargin{:}), '\n'];
out = @(varargin) fprintf(str,varargin{:});

% Optimize for sparsity
function X = optimizeStorage(X,tau)
if nnz(X) < (tau * numel(X))
    if ~issparse(X)
        % Convert to sparse matrix
        X = sparse(X);
    end
else
    if issparse(X)
        % Convert to full matrix
        X = full(X);
    end
end

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
function [D0, B0, type, dr, ddim, fixedD, nIters, order, tau, flag] = ...
                                                        parseInputs(Y,opts)
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

% Standard args
D0     = parseField(opts,'D0',nan);
B0     = parseField(opts,'B0',nan);
type   = parseField(opts,'type','hard');
dr     = parseField(opts,'dr',nan);
ddim   = parseField(opts,'ddim',nan);
fixedD = parseField(opts,'fixedD',false);
nIters = parseField(opts,'nIters',50);
order  = parseField(opts,'order',nan);
tau    = parseField(opts,'tau',-1);
flag   = parseField(opts,'flag',1);

% Expensive defaults
if isnan(D0),    D0    = dctmtx(size(Y,1))';            end
if isnan(B0),    B0    = zeros(size(D0,2),size(Y,2));   end
if isnan(order), order = 1:size(D0,2);                  end

% Parse struct field
function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end

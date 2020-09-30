function [D, Bt, params, stats] = onlineSoupDil(Yt,mu,gamma,varargin)
%
% Syntax:       [D, Bt, params] = onlineSoupDil(Yt,mu,gamma);
%               [D, Bt, params] = onlineSoupDil(Yt,mu,gamma,opts);
%               [D, Bt, params, stats] = onlineSoupDil(Yt,mu,gamma);
%               [D, Bt, params, stats] = onlineSoupDil(Yt,mu,gamma,opts);
%               
% Inputs:       Yt is the p x n data matrix for time t
%               
%               mu >= 0 is the sparsity regularization parameter
%               
%               gamma in [0, 1] is the forgetting factor, where gamma = 0
%               is memoryless, and gamma = 1 is infinite memory
%               
%               [OPTIONAL] opts is a struct containing one or more of the
%               following fields. The default values are in ()
%                   
%                   opts.D0 (dctmtx(p)') is a p x m matrix containing the
%                   initial dictionary for time t (usually the dictionary
%                   from time t - 1). Overcomplete dictionaries (m > p) are
%                   allowed
%                   
%                   opts.B0 (zeros(m,n)) is an m x n matrix containing the
%                   initial sparse codes for time t
%                   
%                   opts.params ([]) is a struct containing the cumulative
%                   parameters returned by the previous call to this 
%                   function. At time t = 1, omit this field or pass
%                   params = [] and it will be initialized appropriately
%                   
%                   opts.type ('hard') can be 'hard' or 'soft' and
%                   specifies whether apply hard or soft thresholding to Bt
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
%                   in which the columns of D and rows of Bt are updated. 
%                   Set order = -1 to randomize the order each iteration
%                   
%                   opts.tau (-1) is the sparsity threshold to use to
%                   optimize the storage of Bt. By default, Bt is a full
%                   matrix, but, when tau > 0, Bt is converted to a sparse
%                   matrix for iterations where nnz(Bt) / numel(Bt) > tau.
%                   Empirically, a good value to try is tau = 0.015
%                   
%                   opts.flag (1) determines what status updates to print
%                   to the command window. The choices are
%                   
%                       flag = 0: no printing
%                       flag = 1: print iteration status
%               
% Outputs:      D is a p x m matrix containing the dictionary at time t
%               
%               Bt is an m x n matrix of sparse codes for time t
%               
%               params is a struct of cumulative parameters for time t (to
%               be passed back into this function at time t + 1)
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
%                   \|Yt - D_k Bt_k\|_F / \|Yt\|_F
%                   
%                   stats.deltaD is a 1 x nIters vector containing the
%                   relative convergence of D at each iteration, defined as
%                   \|D_k - D_{k - 1}\|_F / \|D_{k - 1}\|_F
%                   
%                   stats.deltaB is a 1 x nIters vector containing the
%                   relative convergence of Bt at each iteration, defined:
%                   \|Bt_k - Bt_{k - 1}\|_F / \|Bt_{k - 1}\|_F
%                   
%                   stats.sparsity is a 1 x nIters vector containing the
%                   sparsity, in percent, of Bt at each iteration
%                   
%                   stats.time is a 1 x nIters vector containing the time,
%                   in seconds, required to perform each outer iteration
%               
% Description:  Performs online dictionary learning via the Online Sum of
%               OUter Products DIctionary Learning (Online SOUP-DIL)
%               problem. That is, when type == 'hard':
%               
%               min_{D,Bt} 0.5 \sum_{j=1}^t gamma^{t-j} \|Yj - D Bj\|_F^2 +
%                          0.5 \mu^2 \|Bt\|_0
%               
%               subject to \|D(:,k)\|_2 = 1
%                          rank(R(D(:,k)) <= dr
%               
%               Or, when type == 'soft':
%               
%               min_{D,Bt} 0.5 \sum_{j=1}^t gamma^{t-j} \|Yj - D Bj\|_F^2 +
%                          \mu \|Bt\|_1
%               
%               subject to \|D(:,k)\|_2 = 1
%                          rank(R(D(:,k)) <= dr
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         January 25, 2017
%               May 15, 2017
%

% Parse inputs
[D, Bt, t, Et, Ft, Gt, type, dr, ddim, fixedD, nIters, order, tau, ...
                                       flag] = parseInputs(Yt,varargin{:});
RANDOMIZE_ORDER  = isequal(order,-1);
PRINT_STATS      = (flag > 0);
COMPUTE_STATS    = PRINT_STATS || (nargout == 4);
UPDATE_PARAMS    = (nargout >= 3);
USE_HARD_THRESH  = strcmpi(type,'hard');
OPTIMIZE_STORAGE = (tau > 0);
LOW_RANK_ATOMS   = (dr < min(ddim));
OPTIMIZE_DICT    = ~fixedD;

% Initialize stats
if COMPUTE_STATS
    % Cost function
    Rtfit = @(D) 0.5 * gamma * real(trace(Et) + ...
                                    trace(D' * (D * Gt - 2 * Ft)));
    Ytfit = @(D,Bt) 0.5 * norm(vec(Yt - D * Bt))^2;
    if strcmpi(type,'hard')
        % Ell-0 regularization
        Btfit = @(Bt) 0.5 * mu^2 * nnz(Bt);
    elseif strcmpi(type,'soft')
        % Ell-1 regularization
        Btfit = @(Bt) mu * norm(vec(Bt),1);
    else
        % Unsupported shrinkage
        error('Unsupported shrinkage type ''%s''',type);
    end
    Psi = @(D,Bt) Rtfit(D) + Ytfit(D,Bt) + Btfit(Bt);
    %Psi = @(D,Bt) Ytfit(D,Bt) + Btfit(Bt); % time-t cost only
    
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

% Online SOUP-DIL
if PRINT_STATS
    fprintf('***** Online SOUP-DIL *****\n');
end
[p, m] = size(D);
d0 = [1; zeros(p - 1,1)]; % Zero-reset atom
for it = 1:nIters
    % Initialize iteration
    itimer = tic;
    if RANDOMIZE_ORDER
        % Randomize order
        order = randperm(m);
    end
    if OPTIMIZE_STORAGE
        % Optimize for sparsity
        Bt = optimizeStorage(Bt,tau);
    end
    if COMPUTE_STATS
        % Save last iterates
        Btlast = Bt;
        Dlast = D;
    end
    
    % Block coordinate descent
    DtYt = D' * Yt; % Precompute for efficiency
    for k = order
        % Update Bt(k,:)
        btk = DtYt(k,:) - (D(:,k)' * D) * Bt + Bt(k,:);
        if USE_HARD_THRESH
            % Hard thresholding
            Mk = (abs(btk) > mu);
            btk(~Mk) = 0;
        else
            % Soft thresholding
            btk = sign(btk) .* max(abs(btk) - mu,0);
            Mk = (btk ~= 0);            
        end
        Bt(k,:) = btk;
        ANY_MK = any(Mk);
        
        % Update D(:,k)
        if OPTIMIZE_DICT
            if ANY_MK
                btMk = btk(Mk)';
                ftk = gamma * Ft(:,k) + Yt(:,Mk) * btMk;
                gtk = gamma * Gt(:,k) + Bt(:,Mk) * btMk;                
                dk = ftk - D * gtk + gtk(k) * D(:,k);
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
        cost(it)     = Psi(D,Bt);
        repError(it) = computeNRMSE(D * Bt,Yt);
        deltaD(it)   = computeNRMSE(D,Dlast);
        deltaB(it)   = computeNRMSE(Bt,Btlast);
        sparsity(it) = 100 * (nnz(Bt) / numel(Bt));
        time(it)     = toc(itimer);
        if PRINT_STATS
            out(it,cost(it),repError(it),deltaD(it),deltaB(it), ...
                sparsity(it),time(it));
        end
    end
end

% Update cumulative parameters
if UPDATE_PARAMS
    params.t  = t + 1;
    params.Et = gamma * Et + Yt * Yt';
    params.Ft = gamma * Ft + Yt * Bt';
    params.Gt = gamma * Gt + Bt * Bt';
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
function [D0, B0, t, Et, Ft, Gt, type, dr, ddim, fixedD, nIters, order, ...
                                          tau, flag] = parseInputs(Yt,opts)
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

% Standard args
D0     = parseField(opts,'D0',nan);
B0     = parseField(opts,'B0',nan);
params = parseField(opts,'params',[]);
type   = parseField(opts,'type','hard');
dr     = parseField(opts,'dr',nan);
ddim   = parseField(opts,'ddim',nan);
fixedD = parseField(opts,'fixedD',false);
nIters = parseField(opts,'nIters',50);
order  = parseField(opts,'order',nan);
tau    = parseField(opts,'tau',-1);
flag   = parseField(opts,'flag',1);

% Expensive defaults
if isnan(D0),    D0    = dctmtx(size(Yt,1))';           end
if isnan(B0),    B0    = zeros(size(D0,2),size(Yt,2));  end
if isnan(order), order = 1:size(D0,2);                  end

% Parse params
if isstruct(params)
    t  = params.t;
    Et = params.Et;
    Ft = params.Ft;
    Gt = params.Gt;
else
    t  = 0;
    Et = zeros(size(D0,1));
    Ft = zeros(size(D0));
    Gt = zeros(size(D0,2));
end

% Parse struct field
function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end

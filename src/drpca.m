function [L, S, D, B, stats] = drpca(Y,lambdaL,lambdaS,lambdaB,varargin)
%
% Syntax:       [L, S, D, B] = drpca(Y,lambdaL,lambdaS,lambdaB);
%               [L, S, D, B] = drpca(Y,lambdaL,lambdaS,lambdaB,opts);
%               [L, S, D, B, stats] = drpca(Y,lambdaL,lambdaS,lambdaB);
%               [L, S, D, B, stats] = drpca(Y,lambdaL,lambdaS,lambdaB,opts);
%               
% Inputs:       Y is an m x n data matrix
%               
%               lambdaL >= 0 is the low-rank regularization parameter.
%               Note: lambdaL is ignored when ~isnan(r)
%               
%               lambdaS >= 0 is the patch sparsity regularization parameter
%               
%               lambdaB >= 0 is the sparse coding regularization parameter.
%               Alternatively, lambdaB can be a 1 x nIters vector
%               specifying a different mu value for each outer iteration
%               
%               [OPTIONAL] opts is a struct containing one or more of the
%               following fields. The default values are in ()
%                   
%                   opts.A (1) is an m x p system matrix
%                   
%                   opts.M (1) is an m x n weight matrix
%                   
%                   opts.p (1) is the ell-p shrinkage parameter to use for
%                   the L updates. p = 1 corresponds to SVT, and p = 0
%                   corresponds to a rank penalty
%                   
%                   opts.r (nan) is the OptShrink rank parameter. When a
%                   non-nan r is specified, lambdaL and p are ignored and
%                   OptShrink is used in place of SVT/ell-p shrinkage for
%                   all L updates
%                   
%                   opts.sdim (size(S0)) are the dimensions of S to use
%                   when extracting patches. Note that sdim need not be the
%                   same as size(S) as long as the two representations are
%                   are interchangable via reshape(). For example, you
%                   could specify sdim = [p, n1, n2] where n1 * n2 = n
%                   
%                   opts.pdim ([8, 8]) is a vector with numel(sdim)
%                   elements specifying the patch sizes to extract from S
%                   
%                   opts.pgap ([2, 2]) is a vector with numel(sdim) 
%                   elements specfiying the patch strides (i.e., shifts) to
%                   use along each dimension of S
%                   
%                   opts.np (0) is the number of (randomly extracted)
%                   patches of L0 + S0 to add to the initial dictionary, D0
%                   
%                   opts.type ('hard') can be {'hard','soft'} and specifies
%                   whether to use hard or soft thresholding to regularize
%                   the sparse codes B
%                   
%                   opts.dr (nan) is a rank constraint on the dictionary
%                   atoms. Note: both dr and ddim must be non-nan to apply
%                   a rank constraint
%                   
%                   opts.ddim ([prod(pdim(1:(end - 1))), pdim(end)]) is a
%                   1 x 2 vector describing how to reshape the dictionary
%                   atoms into a matrix before applying the rank constraint
%                   Note: both dr and ddim must be non-nan to apply a rank
%                   constraint
%                   
%                   opts.unitaryD (false) determines whether to enforce a
%                   unitary constraint on the dictionary (true) or not
%                   (false). Note: if unitaryD == true, the np, dr, and
%                   ddim arguments are ignored
%                   
%                   opts.fixedD (false) determines whether to fix the
%                   initial dictionary (true) or learn it (false)
%                   
%                   opts.nIters (50) is the number of outer iterations to
%                   perform
%                   
%                   opts.nItersDB (1) is the number of inner (D,B) updates
%                   to perform per outer iteration
%                   
%                   opts.nItersLS (5) is the number of inner (L,S) updates
%                   to perform per outer iteration
%                   
%                   opts.L0 (A' * Y) is a p x n matrix specifying the
%                   initial L iterate
%                   
%                   opts.S0 (zeros(p,n)) is a p x n matrix specifying the
%                   initial S iterate
%                   
%                   opts.D0 (dctmtx(b)') is a b x c matrix specifying the
%                   initial dictionary. Undercomplete (b > c) and
%                   overcomplete (b < c) dictionaries are allowed, unless
%                   unitaryD == true, in which case we must have b = c
%                   
%                   opts.B0 (zeros(c,d)) is an c x d matrix specifying the
%                   initial sparse codes
%                   
%                   opts.Xtrue (nan) is the ground truth X = L + S matrix
%                   to use for NRMSE calculations
%                   
%                   opts.NRMSEfcn (@computeNRMSE) is the function to use
%                   when computing the NRMSE of X = L + S after each
%                   iteration
%                   
%                   opts.accel (true) specifies whether to use Nesterov's
%                   acceleration scheme for the (L,S) updates
%                   
%                   opts.tau 0.99 / ((1 + accel) * norm(A)^2) is the step
%                   size for the (L,S) updates, and should satisfy
%                   
%                       tau <= 1 / (2 * norm(A)^2), when accel = true
%                       tau <  1 /      norm(A)^2 , when accel = false
%                   
%                   opts.flag (1) determines what status updates to print
%                   to the command window. The choices are
%                   
%                       flag = 0: no printing
%                       flag = 1: print outer iteration status
%                       flag = 2: print inner and outer iteration status
%               
% Outputs:      L is a p x n matrix containing the estimated low-rank
%               component
%               
%               S is a p x n matrix containing the estimated patch-sparse
%               component
%               
%               D is a b x c matrix, where b = prod(pdim), c = size(D0,2), 
%               whose columns contain the vectorized dictionary elements
%               
%               B is a c x d matrix, where d = total # patches, whose
%               columns contain the sparse codes for the patches of S
%               
%               stats is a struct containing the following fields:
%               
%                   stats.nIters is the number of outer iterations
%                   performed
%                   
%                   stats.cost is a 1 x nIters vector containing the cost
%                   function values at each iteration
%                   
%                   stats.nrmse is a 1 x nIters vector containing the NRMSE
%                   of X = L + S with respect to Xtrue at each iteration
%                   
%                   stats.delta is a 1 x nIters vector containing the
%                   relative convergence of X = L + S at each iteration,
%                   defined as \|X_k - X_{k - 1}\|_F / \|X_{k - 1}\|_F
%                   
%                   stats.sparsity is a 1 x nIters vector containing the
%                   sparsity, in percent, of B at each iteration
%                   
%                   stats.time is a 1 x nIters vector containing the time,
%                   in seconds, required to perform each outer iteration
%               
% Description:  Solves one of the following Dictionary Robust PCA (DRPCA)
%               problems:
%               
%               When type == 'hard':
%               
%             \argmin_{L,S,D,B} 0.5\|\sqrt{M} \odot (Y - A(L + S))\|_F^2 +
%                                  \lambda_L\|L\|_{\star}  +
%                 \lambda_S (0.5\|P(S) - DB\|_F^2 + 0.5\lambda_B^2\|B\|_0)
%               
%               When type == 'soft':
%               
%             \argmin_{L,S,D,B} 0.5\|\sqrt{M} \odot (Y - A(L + S))\|_F^2 +
%                                  \lambda_L\|L\|_{\star}  +
%                 \lambda_S (0.5\|P(S) - DB\|_F^2 + \lambda_B\|B\|_1)
%               
%               subject to
%                   \|D(:,k)\|_2 = 1             (unitaryD == false)
%                   rank(R(D(:,k)) <= dr         (~isnan(dr))
%               or
%                   D' D = I                     (unitaryD == true)
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         June 17, 2016
%               November 6, 2016
%               May 15, 2017
%               June 6, 2017
%

% Parse inputs
[A, M, p, r, sdim, pdim, pgap, np, type, dr, ddim, unitaryD, fixedD, ...
 nIters, nItersDB, nItersLS, L, S, D, B, Xtrue, NRMSEfcn, accel, tau, ...
                                        flag] = parseInputs(Y,varargin{:});
if isscalar(lambdaB)
    lambdaB = repmat(lambdaB,1,nIters);
end
PRINT_STATS   = (flag > 0);
COMPUTE_STATS = PRINT_STATS || (nargout == 5);
USE_OPTSHRINK = ~isnan(r);
UNITARY_DICT  = unitaryD;

% Initialize stats
if COMPUTE_STATS
    % Cost function
    Xfit  = @(X) 0.5 * norm(vec(sqrt(M) .* (Y - A * X)))^2;
    if p == 0
        % Rank penalty
        Lfit = @(sL) 0.5 * lambdaL^2 * nnz(sL);
    else
        % Ell-p spectral regularization
        Lfit = @(sL) lambdaL * sum(sL.^p);
    end
    Pfit = @(P,DB,B) 0.5 * lambdaS * norm(vec(P - DB))^2;
    if strcmpi(type,'hard')
        % Ell-0 regularization
        Bfit = @(B) 0.5 * lambdaS * lambdaB(end)^2 * nnz(B);
    elseif strcmpi(type,'soft')
        % Ell-1 regularization
        Bfit = @(B) lambdaS * lambdaB(end) * norm(vec(B),1);        
    else
        % Unsupported shrinkage
        error('Unsupported shrinkage type ''%s''',type);
    end
    Psi = @(X,sL,P,DB,B) Xfit(X) + Lfit(sL) + Pfit(P,DB,B) + Bfit(B);
    
    % Stats-printing function
    iterFmt = sprintf('%%0%dd',ceil(log10(nIters + 1)));
    out     = printFcn('Iter'    ,iterFmt ,'cost' ,'%.2f', ...
                       'nrmse'   ,'%.3f'  ,'delta','%.3e', ...
                       'sparsity','%.2f%%','time' ,'%.2fs');
    
    % Initialize stats
    X        = L + S;
    cost     = nan(1,nIters);
    nrmse    = nan(1,nIters);
    delta    = nan(1,nIters);
    sparsity = nan(1,nIters);
    time     = nan(1,nIters);
end

% Initialize patches
[Mp, ~, Np] = patchInds(sdim,pdim,pgap); % Patch indices
Np = reshape(Np,size(S));                % Patch counts
P  = S(Mp);                              % Initial patches

% Modify dictionary, if requested
if ~UNITARY_DICT
    if np > 0
        % Append random patches
        mp   = size(Mp,2);
        idx  = randperm(mp);
        inds = Mp(:,idx(1:min(np,mp)));
        D    = [D, (L(inds) + S(inds))];
    elseif np < 0
        % Remove random patches
        mp  = size(D,2);
        idx = randperm(mp);
        D   = D(:,idx(1:max(1,mp + np)));
    end
end

% Initialize dictionary learning opts
optsDB.type   = type;
optsDB.fixedD = fixedD;
optsDB.nIters = nItersDB;
optsDB.flag   = flag - 1;
if ~UNITARY_DICT
    % soupDil-specific opts
    optsDB.dr   = dr;
    optsDB.ddim = ddim;
end

% Initialize patchRPCA opts
optsLS.A        = A;
optsLS.M        = M;
optsLS.p        = p;
optsLS.r        = r;
optsLS.nIters   = nItersLS;
optsLS.Xtrue    = Xtrue;
optsLS.NRMSEfcn = NRMSEfcn;
optsLS.accel    = accel;
optsLS.tau      = tau;
optsLS.flag     = flag - 1;

% Dictionary Robust PCA
if PRINT_STATS
    fprintf('***** Dictionary Robust PCA *****\n');
end
if isnan(B)
    % Initialize sparse codes
    B = zeros(size(D,2),size(P,2));
end
for it = 1:nIters
    % Initialize iteration
    itimer = tic;
    
    % (D,B) update(s)
    optsDB.D0 = D;
    optsDB.B0 = B;
    if UNITARY_DICT
        % unitaryDil
        [D, B] = unitaryDil(P,lambdaB(it),optsDB);
    else
        % soupDil
        [D, B] = soupDil(P,lambdaB(it),optsDB);
    end
    DB = D * B;
    
    % (L,S) update(s)
    optsLS.L0 = L;
    optsLS.S0 = S;
    [L, S, sL] = patchRPCA(Y,DB,lambdaL,lambdaS,Mp,Np,optsLS);
    P = S(Mp);
    
    % Record stats
    if COMPUTE_STATS
        Xlast = X;
        X     = L + S;
        if ~USE_OPTSHRINK
            cost(it) = Psi(X,sL,P,DB,B);
        end
        nrmse(it)    = NRMSEfcn(X,Xtrue);
        delta(it)    = computeNRMSE(X,Xlast);
        sparsity(it) = 100 * (nnz(B) / numel(B));
        time(it)     = toc(itimer);
        if PRINT_STATS
            out(it,cost(it),nrmse(it),delta(it),sparsity(it),time(it)); 
        end
    end
end

% Return stats
if COMPUTE_STATS
    stats.nIters   = nIters;
    stats.nrmse    = nrmse;
    stats.cost     = cost;
    stats.delta    = delta;
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
function [A, M, p, r, sdim, pdim, pgap, np, type, dr, ddim, unitaryD, ...
          fixedD, nIters, nItersDB, nItersLS, L0, S0, D0, B0, Xtrue, ...
          NRMSEfcn, accel, tau, flag] = parseInputs(Y,opts)
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

% Standard args
A        = parseField(opts,'A',1);
M        = parseField(opts,'M',1);
p        = parseField(opts,'p',1);
r        = parseField(opts,'r',nan);
%sdim dependent
pdim     = parseField(opts,'pdim',[8, 8]);
pgap     = parseField(opts,'pgap',[2, 2]);
np       = parseField(opts,'np',0);
type     = parseField(opts,'type','hard');
dr       = parseField(opts,'dr',nan);
%ddim dependent
unitaryD = parseField(opts,'unitaryD',false);
fixedD   = parseField(opts,'fixedD',false);
nIters   = parseField(opts,'nIters',50);
nItersDB = parseField(opts,'nItersDB',1);
nItersLS = parseField(opts,'nItersLS',5);
L0       = parseField(opts,'L0',nan);
S0       = parseField(opts,'S0',nan);
D0       = parseField(opts,'D0',nan);
B0       = parseField(opts,'B0',nan);
Xtrue    = parseField(opts,'Xtrue',nan);
NRMSEfcn = parseField(opts,'NRMSEfcn',@computeNRMSE);
accel    = parseField(opts,'accel',true);
tau      = parseField(opts,'tau',nan);
flag     = parseField(opts,'flag',1);

% Expensive defaults
if isnan(L0),  L0  = A' * Y;                        end
if isnan(S0),  S0  = zeros(size(L0));               end
if isnan(D0),  D0  = dctmtx(prod(pdim))';           end
if isnan(tau), tau = 1 / ((1 + accel) * norm(A)^2); end

% Dependent args
sdim = parseField(opts,'sdim',size(S0));
ddim = parseField(opts,'ddim',[prod(pdim(1:(end - 1))), pdim(end)]);

% Parse struct field
function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end

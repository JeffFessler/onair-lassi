function [Xt, D, Bt, params, stats] = onlineDls(Yt,lambda,mu,gamma,varargin)
%
% Syntax:       [Xt, D, Bt, params] = onlineDls(Yt,lambda,mu,gamma);
%               [Xt, D, Bt, params] = onlineDls(Yt,lambda,mu,gamma,opts);
%               [Xt, D, Bt, params, stats] = onlineDls(Yt,lambda,mu,gamma);
%               [Xt, D, Bt, params, stats] = onlineDls(Yt,lambda,mu,gamma,opts);
%               
% Inputs:       Yt is an m x n data matrix
%               
%               lambda >= 0 is the dictionary regularization parameter
%               
%               mu >= 0 is the sparsity regularization parameter.
%               Alternatively, mu can be a 1 x nIters vector specifying a
%               different mu value for each outer iteration
%               
%               gamma in [0, 1] is the forgetting factor, where gamma = 0
%               is memoryless, and gamma = 1 is infinite memory
%               
%               [OPTIONAL] opts is a struct containing one or more of the
%               following fields. The default values are in ()
%                   
%                   opts.A (1) is an m x p system matrix
%                   
%                   opts.M (1) is an m x n weight matrix
%                   
%                   opts.xdim (size(X0)) are the dimensions of Xt to use
%                   when extracting patches. Note that xdim need not be the
%                   same as size(Xt) as long as the two representations are
%                   are interchangable via reshape(). For example, you
%                   could specify xdim = [p, n1, n2] where n1 * n2 = n
%                   
%                   opts.pdim ([8, 8]) is a vector with numel(xdim)
%                   elements specifying the patch sizes to extract from Xt
%                   
%                   opts.pgap ([2, 2]) is a vector with numel(xdim) 
%                   elements specfiying the patch strides (i.e., shifts) to
%                   use along each dimension of Xt
%                   
%                   opts.type ('hard') can be {'hard','soft'} and specifies
%                   whether to use hard or soft thresholding to regularize
%                   the sparse codes Bt
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
%                   (false). Note: if unitaryD == true, the dr and ddim
%                   arguments are ignored
%                   
%                   opts.fixedD (false) determines whether to fix the
%                   initial dictionary (true) or learn it (false)
%                   
%                   opts.nIters (50) is the number of outer iterations to
%                   perform
%                   
%                   opts.nItersDB (1) is the number of inner (D,Bt) updates
%                   to perform per outer iteration
%                   
%                   opts.nItersX (5) is the number of inner Xt updates to
%                   perform per outer iteration. Note that, when A = 1, the
%                   exact Xt update is always computed and nItersX is
%                   ignored
%                   
%                   opts.X0 (A' * Yt) is a p x n matrix specifying the
%                   initial Xt iterate
%                   
%                   opts.D0 (dctmtx(b)') is a b x c matrix specifying the
%                   initial dictionary for time t (usually, the dictionary
%                   from time t - 1). Undercomplete (b > c) and
%                   overcomplete (b < c) dictionaries are allowed, unless
%                   unitaryD == true, in which case we must have b = c
%                   
%                   opts.B0 (zeros(c,d)) is an c x d matrix specifying 
%                   the initial sparse codes for time t
%                   
%                   opts.params ([]) is a struct containing the cumulative
%                   parameters returned by the previous call to this 
%                   function. At time t = 1, omit this field or pass
%                   params = [] and it will be initialized appropriately
%                   
%                   opts.Xtrue (nan) is the ground truth Xt matrix to use 
%                   for NRMSE calculations
%                   
%                   opts.NRMSEfcn (@computeNRMSE) is the function to use
%                   when computing the NRMSE of Xt after each iteration
%                   
%                   opts.accel (true) specifies whether to use Nesterov's
%                   acceleration scheme for the Xt updates
%                   
%                   opts.tau ((0.99 + ~accel) / norm(A)^2) is the step size
%                   parameter for the Xt updates, and should satisfy
%                   
%                       tau <= 1 / norm(A)^2, when accel = true
%                       tau <  2 / norm(A)^2, when accel = false
%                   
%                   opts.flag (1) determines what status updates to print
%                   to the command window. The choices are
%                   
%                       flag = 0: no printing
%                       flag = 1: print outer iteration status
%                       flag = 2: print inner and outer iteration status
%               
% Outputs:      Xt is the estimated p x n matrix for time t
%               
%               D is a b x c matrix, where b = prod(pdim), c = size(D0,2), 
%               whose columns contain the dictionary atoms for time t
%               
%               Bt is a c x d matrix, where d = total # patches, whose
%               columns are sparse codes for the patches of Xt at time t
%               
%               params is a struct of cumulative parameters for time t (to
%               be passed back into this function at time t + 1)
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
%                   of Xt with respect to Xtrue at each iteration
%                   
%                   stats.deltaX is a 1 x nIters vector containing the
%                   relative convergence of Xt at each iteration, defined:
%                   \|Xt_{k + 1} - Xt_k\|_F / \|Xt_k\|_F
%                   
%                   stats.deltaD is a 1 x nIters vector containing the
%                   relative convergence of D at each iteration, defined as
%                   \|D_{k + 1} - D_k\|_F / \|D_k\|_F
%                   
%                   stats.deltaB is a 1 x nIters vector containing the
%                   relative convergence of Bt at each iteration, defined:
%                   \|Bt_{k + 1} - Bt_k\|_F / \|Bt_k\|_F
%                   
%                   stats.sparsity is a 1 x nIters vector containing the
%                   sparsity, in percent, of Bt at each iteration
%                   
%                   stats.time is a 1 x nIters vector containing the time,
%                   in seconds, required to perform each outer iteration
%               
% Description:  Solves one of the following Online Dictionary Least-Sqaures
%               (Online DLS) problems:
%               
%               When type == 'hard':
%               
%       \min_{Xt,D,Bt} 0.5 \|\sqrt{M} \odot (Yt - A Xt)\|_F^2 + 
%             \lambda (0.5 \sum_{j=1}^t gamma^{t-j} \|Pj(Xj) - D Bj\|_F^2 +
%                      0.5 \mu^2 \|Bt\|_0)
%               
%               When type == 'soft':
%               
%       \min_{Xt,D,Bt} 0.5 \|\sqrt{M} \odot (Yt - A Xt)\|_F^2 + 
%             \lambda (0.5 \sum_{j=1}^t gamma^{t-j} \|Pj(Xj) - D Bj\|_F^2 +
%                          \mu \|Bt\|_1)
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
% Date:         January 25, 2017
%               May 15, 2017
%               June 6, 2017
%

% Parse inputs
[mu, A, M, xdim, pdim, pgap, type, dr, ddim, unitaryD, fixedD, nIters, ...
 nItersDB, nItersX, Xt, D, Bt, params, Xtrue, NRMSEfcn, accel, tau, ...
                                    flag] = parseInputs(Yt,mu,varargin{:});
PRINT_STATS   = (flag > 0);
UPDATE_PARAMS = (nargout >= 4);
COMPUTE_STATS = PRINT_STATS || (nargout == 5);
UNITARY_DICT  = unitaryD;

% Initialize stats
if COMPUTE_STATS
    % Cost function
    Ytfit = @(Xt) 0.5 * norm(vec(sqrt(M) .* (Yt - A * Xt)))^2;
    Rtfit = @(D) 0.5 * gamma * real(trace(params.Et) + ...
                              trace(D' * (D * params.Gt - 2 * params.Ft)));
    Xtfit = @(Pt,DBt) 0.5 * norm(vec(Pt - DBt))^2;
    if strcmpi(type,'hard')
        % Ell-0 regularization
        Btfit = @(Bt) 0.5 * mu(end)^2 * nnz(Bt);
    elseif strcmpi(type,'soft')
        % Ell-1 regularization
        Btfit = @(Bt) mu(end) * norm(vec(Bt),1);
    else
        % Unsupported shrinkage
        error('Unsupported shrinkage type ''%s''',type);
    end
    Psi = @(Xt,Pt,D,DBt,Bt) Ytfit(Xt) + ...
                           lambda * (Rtfit(D) + Xtfit(Pt,DBt) + Btfit(Bt));
    %{
    % time-t cost only 
    Psi = @(Xt,Pt,DBt,Bt) Ytfit(Xt) + lambda * (Xtfit(Pt,DBt) + Btfit(Bt));
    %}
    
    % Stats-printing function
    iterFmt = sprintf('%%0%dd',ceil(log10(nIters + 1)));
    out     = printFcn('Iter'    ,iterFmt ,'cost'  ,'%.2f', ...
                       'nrmse'   ,'%.3f'  ,'deltaX','%.3e', ...
                       'deltaD'  ,'%.3e'  ,'deltaB','%.3e', ...
                       'sparsity','%.2f%%','time'  ,'%.2fs');
    
    % Initialize stats
    cost     = nan(1,nIters);
    nrmse    = nan(1,nIters);
    deltaX   = nan(1,nIters);
    deltaD   = nan(1,nIters);
    deltaB   = nan(1,nIters);
    sparsity = nan(1,nIters);
    time     = nan(1,nIters);
end

% Initialize patches
[Mp, ~, Np] = patchInds(xdim,pdim,pgap); % Patch indices
Np = reshape(Np,size(Xt));               % Patch counts
Pt = Xt(Mp);                             % Initial patches

% Initialize online dictionary learning opts
optsDB.params = params;
optsDB.type   = type;
optsDB.fixedD = fixedD;
optsDB.nIters = nItersDB;
optsDB.flag   = flag - 1;
if ~UNITARY_DICT
    % onlineSoupDil-specific opts
    optsDB.dr   = dr;
    optsDB.ddim = ddim;
end

% Initialize patchLS opts
optsX.A        = A;
optsX.M        = M;
optsX.nIters   = nItersX;
optsX.Xtrue    = Xtrue;
optsX.NRMSEfcn = NRMSEfcn;
optsX.accel    = accel;
optsX.tau      = tau;
optsX.flag     = flag - 1;

% Online Dictionary Least-Squares
if PRINT_STATS
    fprintf('***** Online Dictionary Least-Squares *****\n');
end
if isnan(Bt)
    % Initialize sparse codes
    Bt = zeros(size(D,2),size(Pt,2));       
end
for it = 1:nIters
    % Initialize iteration
    itimer = tic;
    if COMPUTE_STATS
        % Save last iterate
        Xlast = Xt;
        Dlast = D;
        Blast = Bt;
    end
    
    % (D,Bt) update(s)
    optsDB.D0 = D;
    optsDB.B0 = Bt;
    if UNITARY_DICT
        % onlineUnitaryDil
        [D, Bt] = onlineUnitaryDil(Pt,mu(it),gamma,optsDB);
    else
        % onlineSoupDil
        [D, Bt] = onlineSoupDil(Pt,mu(it),gamma,optsDB);
    end
    DBt = D * Bt;
    
    % Xt update(s)
    optsX.X0 = Xt;
    Xt = patchLS(Yt,DBt,lambda,Mp,Np,optsX);
    Pt = Xt(Mp);
    
    % Record stats
    if COMPUTE_STATS
        cost(it)     = Psi(Xt,Pt,D,DBt,Bt);
        nrmse(it)    = NRMSEfcn(Xt,Xtrue);
        deltaX(it)   = computeNRMSE(Xt,Xlast);
        deltaD(it)   = computeNRMSE(D,Dlast);
        deltaB(it)   = computeNRMSE(Bt,Blast);
        sparsity(it) = 100 * (nnz(Bt) / numel(Bt));
        time(it)     = toc(itimer);
        if PRINT_STATS
            out(it,cost(it),nrmse(it),deltaX(it),deltaD(it),deltaB(it), ...
                sparsity(it),time(it)); 
        end
    end
end

% Update cumulative parameters
if UPDATE_PARAMS
    params.t  = params.t + 1;
    params.Et = gamma * params.Et + Pt * Pt';
    params.Ft = gamma * params.Ft + Pt * Bt';
    params.Gt = gamma * params.Gt + Bt * Bt';
end

% Return stats
if COMPUTE_STATS
    stats.nIters   = nIters;
    stats.cost     = cost;
    stats.nrmse    = nrmse;
    stats.deltaX   = deltaX;
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
function [mu, A, M, xdim, pdim, pgap, type, dr, ddim, unitaryD, fixedD, ...
          nIters, nItersDB, nItersX, X0, D0, B0, params, Xtrue, ...
                      NRMSEfcn, accel, tau, flag] = parseInputs(Yt,mu,opts)
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

% Standard args
A        = parseField(opts,'A',1);
M        = parseField(opts,'M',1);
%xdim dependent
pdim     = parseField(opts,'pdim',[8, 8]);
pgap     = parseField(opts,'pgap',[2, 2]);
type     = parseField(opts,'type','hard');
dr       = parseField(opts,'dr',nan);
%ddim dependent
unitaryD = parseField(opts,'unitaryD',false);
fixedD   = parseField(opts,'fixedD',false);
nIters   = parseField(opts,'nIters',50);
nItersDB = parseField(opts,'nItersDB',1);
nItersX  = parseField(opts,'nItersX',5);
X0       = parseField(opts,'X0',nan);
D0       = parseField(opts,'D0',nan);
B0       = parseField(opts,'B0',nan);
params   = parseField(opts,'params',[]);
Xtrue    = parseField(opts,'Xtrue',nan);
NRMSEfcn = parseField(opts,'NRMSEfcn',@computeNRMSE);
accel    = parseField(opts,'accel',true);
tau      = parseField(opts,'tau',nan);
flag     = parseField(opts,'flag',1);

% Expensive defaults
if isnan(X0),    X0  = A' * Yt;                     end
if isnan(D0),    D0  = dctmtx(prod(pdim))';         end
if isnan(tau),   tau = (0.99 + ~accel) / norm(A)^2; end
if isscalar(mu), mu  = repmat(mu,1,nIters);         end

% Dependent args
xdim = parseField(opts,'xdim',size(X0));
ddim = parseField(opts,'ddim',[prod(pdim(1:(end - 1))), pdim(end)]);

% Parameters
if ~isstruct(params)
    params.t  = 0;
    params.Et = zeros(size(D0,1));
    params.Ft = zeros(size(D0));
    params.Gt = zeros(size(D0,2));
end

% Parse struct field
function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end

function [X, stats] = ktSLR(Y,A,p,mu1,mu2,X0,varargin)
%
% Syntax:       X = ktSLR(Y,A,p,mu1,mu2,X0);
%               X = ktSLR(Y,A,p,mu1,mu2,X0,opts);
%               [X, stats] = ktSLR(Y,A,p,mu1,mu2,X0);
%               [X, stats] = ktSLR(Y,A,p,mu1,mu2,X0,opts);
%               
% Inputs:       Y is a ny x nx x nt x nc matrix containing k-space data
%               
%               A is the system matrix
%               
%               p is the matrix ell-p norm parameter
%               
%               mu1 is the low-rank regularization parameter
%               
%               mu2 is the TV regularization parameter
%               
%               X0 is a ny x nx x nt matrix containing the initial X
%               iterate
%               
%               [OPTIONAL] opts is a struct containing one or more of the
%               following fields. The default values are in ()
%                   
%                   opts.nItersO (5) is the number of outer iterations to
%                   perform
%                   
%                   opts.nItersI (10) is the number of inner iterations to
%                   perform
%                   
%                   opts.nItersCG (5) is the number of outer iterations to
%                   perform
%                   
%                   opts.beta1 (1) is the initial low-rank continuation
%                   parameter
%                   
%                   opts.beta2 (1) is the initial TV continuation
%                   parameter
%                   
%                   opts.beta1rate (50) is the low-rank continuation
%                   multiplicative increment
%                   
%                   opts.beta2rate (50) is the TV continuation
%                   multiplicative increment
%                   
%                   opts.stepSize ([1, 1, 0.303]) are the finite difference
%                   step sizes
%                   
%                   opts.Xtrue (nan) is the ground truth X matrix to use
%                   for NRMSE calculations
%                   
% Outputs:      X is an ny x nx x nt matrix containing the reconstructed
%               data
%               
%               stats is a statistics struct containing the following
%               fields:
%               
%                   stats.nIters is the number of iterations performed
%                   
%                   stats.cost is a 1 x nIters vector containing the value
%                   of the Robust PCA cost function at each iteration
%                   
%                   stats.nrmse is a 1 x nIters vector containing the NRMSE
%                   of X = L + S with respect to Xtrue at each iteration
%                   
%                   stats.delta is a 1 x nIters vector containing the
%                   relative convergence of X = L + S at each iteration,
%                   defined as \|X_{k + 1} - X_k\|_F / \|X_k\|_F
%                   
%                   stats.time is a 1 x nIters vector containing the time,
%                   in seconds, required to perform each iteration
%               
% Description:  Solves the k-t SLR optimization problem
%               
%               X = \argmin_{X} 0.5\|Y - A(X))\|_F^2 +
%                               \mu_1 \|X\|_{p} + \mu_2 TV(X)
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         May 27, 2016
%

% Parse inputs
[nItersO, nItersI, nItersCG, beta1, beta2, beta1rate, beta2rate, ...
                               stepSize, Xtrue] = parseInputs(varargin{:});
[ny, nx, nt]  = size(X0); 
nIters        = nItersO * nItersI;
COMPUTE_STATS = (nargout == 2);

% Generate function handles
Afcn  = @(Z) A * Z;
Atfcn = @(Z) reshape(A' * Z,ny,nx,nt);
[Dfcn, Dtfcn] = defDDt(stepSize);

% Initialize stats
if COMPUTE_STATS
    % Cost function
    Psi = @(X,DX) norm(vec(Afcn(X) - Y))^2 + ...
            mu1 * ellpnorm(reshape(X,[],nt),p) + ...
            mu2 * sum(vec(abs(sqrt(DX{1}.^2 + DX{2}.^2 + DX{3}.^2))));
    
    % Stats-printing function
    iterFmt = sprintf('%%0%dd',ceil(log10(nIters + 1)));
    out     = printFcn('Iter' ,iterFmt,'cost','%.2f','nrmse','%.3f', ...
                       'delta','%.3e' ,'time','%.2fs');
    
    % Initialize stats
    cost   = nan(1,nIters);
    nrmse  = nan(1,nIters);
    delta  = nan(1,nIters);
    time   = nan(1,nIters);
end

% k-t SLR
fprintf('***** k-t SLR *****\n');
X    = X0;
Lam1 = zeros(ny,nx,nt);
Lam2 = zeros(ny,nx,nt);
Lam3 = zeros(ny,nx,nt);
Lam4 = zeros(ny,nx,nt);
DX   = Dfcn(X);
it   = 0;
for itO = 1:nItersO
    for itI = 1:nItersI
        % Initialize iteration
        if COMPUTE_STATS
            it     = it + 1;
            itimer = tic;
            Xlast  = X;
        end
        
        % W update (TV)
        Z1    = DX{1} + Lam1 / beta2;
        Z2    = DX{2} + Lam2 / beta2;
        Z3    = DX{3} + Lam3 / beta2;
        V     = sqrt(abs(Z1).^2 + abs(Z2).^2 + abs(Z3).^2);
        V(~V) = 1;
        V     = max(0,V - 1 / beta2) ./ V;
        Wx    = Z1 .* V;
        Wy    = Z2 .* V;
        Wt    = Z3 .* V;
        
        % Lambda update (matrix ell-p value shrinkage)
        XL4    = reshape(X + Lam4 / beta1,[],nt);
        Lambda = reshape(ellpshrink(XL4,p,beta1),ny,nx,nt);
        
        % X update (CG)
        X  = cg(Y,X,Afcn,Atfcn,Dfcn,Dtfcn,Lambda,Lam1,Lam2,Lam3,Lam4,Wx,Wy,Wt,mu1,mu2,beta1,beta2,nItersCG);
        DX = Dfcn(X);
        
        % Lagrange multiplier updates
        Lam1 = Lam1 - 1.618 * beta2 * (Wx - DX{1});
        Lam2 = Lam2 - 1.618 * beta2 * (Wy - DX{2});
        Lam3 = Lam3 - 1.618 * beta2 * (Wt - DX{3});
        Lam4 = Lam4 - 1.618 * beta1 * (Lambda - X);
        
        % Record stats
        if COMPUTE_STATS
            cost(it)  = Psi(X,DX);
            nrmse(it) = computeNRMSE(X,Xtrue);        
            delta(it) = computeNRMSE(X,Xlast);
            time(it)  = toc(itimer);
            out(it,cost(it),nrmse(it),delta(it),time(it)); 
        end
        
        %{
        % Check for convergence
        if (itI > 1) && (abs(cost(it) - cost(it - 1)) < 1e-3 * abs(cost(it)))
            disp('Converged!');
            break;
        end
        %}
    end
    
    % Beta updates
    beta1 = beta1 * beta1rate;
    beta2 = beta2 * beta2rate;
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

% Matrix ell-p norm
function np = ellpnorm(X,p)
[~, Su, ~] = givefastSVD(X);
np         = sum(diag(Su).^p) / p;

% Matrix ell-p shrinkage
function Xp = ellpshrink(X,p,beta)
[U, S, V] = givefastSVD(X);
s         = diag(S);
Xp        = U * diag(max(0,s - (s.^(p - 1)) / beta)) * V';

% Vectorize data
function x = vec(X)
x = X(:);

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

% Parse inputs
function [nItersO, nItersI, nItersCG, beta1, beta2, beta1rate, ...
                            beta2rate, stepSize, Xtrue] = parseInputs(opts)
if ~exist('opts','var')
    opts = struct();
end

% Standard args
nItersO   = parseField(opts,'nItersO',5);
nItersI   = parseField(opts,'nItersI',10);
nItersCG  = parseField(opts,'nItersCG',5);
beta1     = parseField(opts,'beta1',1);
beta2     = parseField(opts,'beta2',1);
beta1rate = parseField(opts,'beta1rate',50);
beta2rate = parseField(opts,'beta2rate',50);
stepSize  = parseField(opts,'stepSize',[1, 1, 0.303]);
Xtrue     = parseField(opts,'Xtrue',nan);

% Parse struct field
function val = parseField(stats,field,default)
if isfield(stats,field)
    val = stats.(field);
else
    val = default;
end

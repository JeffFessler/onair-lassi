function LASSI_full_run(parFcn,idx)
% Syntax: LASSI_full_run();
%         LASSI_full_run(parFcn,idx);

% Get parameters
[vars, opts, vars0, opts0] = parFcn();

% Get parameters for this iteration
nP = numel(vars.ps);
nS = numel(vars.seed);
nT = nP * nS;
if ~exist('idx','var') || isempty(idx)
    % Return job count
    fprintf('Total jobs: %d\n',nT);
    return;
elseif (idx < 1) || (idx > nT)
    % Invalid index
    fprintf('Index %d is out of range [1,%d]\n',idx,nT);
    return;
else
    % Valid index
    fprintf('Running index %d/%d\n',idx,nT);
end
[ii, jj] = ind2sub([nP, nS],idx);
ps   = vars.ps(ii);
seed = vars.seed(jj);

% Parse path
[path, name, ext] = fileparts(vars.rawpath);
if ~exist(path,'dir')
    % Make results directory
    mkdir(path);
end
out = sprintf('%s/%s%d%s',path,name,idx,ext);
if exist(out,'file')
    % Output data already exists
    fprintf('Output file ''%s'' already exists\n',out);
    return;
end

% Add dependencies to path
addpath('./deps_lassi');

% Generate undersampled data
[Y, A, ~, Xtrue, Xfft] = generateCardiacPerfData(ps,inf,seed,vars.inpath);
[ny, nx, nt] = size(Xtrue);

% Initialization
if isempty(vars0)
    % FFT-based initialization
    L0 = zeros(ny,nx,nt);
    S0 = Xfft;
elseif isfield(vars0,'lambdaL')
    % Add dependencies to path
    addpath('./deps_rpca');
    
    % Run RPCA
    opts0.A      = A;
    opts0.T      = TempFFT(3,[ny, nx, nt]);
    opts0.r      = vars0.r;
    opts0.nIters = vars0.nIters;
    opts0.L0     = reshape(Xfft,[],nt);
    opts0.S0     = zeros(ny * nx,nt);
    [L0, S0]     = robustPCA(Y,vars0.lambdaL,vars0.lambdaS,opts0);
    if isfield(vars0,'L0SX') && (vars0.L0SX == true)
        % L = 0, S = X
        X0 = L0 + S0;
        L0 = zeros(size(X0));
        S0 = X0;
    end
elseif isfield(vars0,'mu1')
    % Add dependencies to path
    addpath('./deps_ktslr');
    
    % Run k-t SLR
    opts0.nItersO  = vars0.nItersO;
    opts0.nItersI  = vars0.nItersI;
    opts0.nItersCG = vars0.nItersCG;
    x0 = ktSLR(Y,A,vars0.p,vars0.mu1,vars0.mu2,Xfft,opts0);
    
    % Construct LASSI initialization
    L0 = zeros(ny,nx,nt);
    S0 = x0;
else
    % Unrecognized algorithm
    error('Unrecognized algorithm');
end

% Run LASSI
opts.A       = A;
if isfield(vars,'p')
    opts.p   = vars.p;
end
opts.r       = vars.r;
opts.sdim    = [ny, nx, nt];
opts.dr      = vars.dr;
opts.nIters  = vars.nIters;
opts.L0      = reshape(L0,[],nt);
opts.S0      = reshape(S0,[],nt);
opts.Xtrue   = reshape(Xtrue,[],nt);
lB = logspace(log10(vars.lambdaB0),log10(vars.lambdaB),vars.nIters);
[Lhat, Shat, ~, ~, stats] = drpca(Y,vars.lambdaL,vars.lambdaS,lB,opts);%#ok

% Aggregate params
params.ps   = ps;
params.seed = seed; %#ok

% Save results
save(out,'Lhat','Shat','params','stats');
fprintf('DONE\n');

function invivo_LASSI_full_run(parFcn,idx)
% Syntax: invivo_LASSI_full_run();
%         invivo_LASSI_full_run(parFcn,idx);

% Get parameters
[vars, opts, vars0, opts0] = parFcn();

% Get parameters for this iteration
nL = numel(vars.nLines);
nS = numel(vars.seed);
nT = nL * nS;
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
[ii, jj] = ind2sub([nL, nS],idx);
nLines = vars.nLines(ii);
seed   = vars.seed(jj);

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
[Y, A, ~, Xtrue, Xfft] = generateInvivoData(nLines,inf,seed,vars.inpath);
[ny, nx, nt] = size(Xtrue);

% Initialization
if isfield(vars0,'lambdaL')
    % Add dependencies to path
    addpath('./deps_rpca');
    
    % Run RPCA
    opts0.A      = A;
    opts0.T      = TempFFT(3,[ny, nx, nt]);
    opts0.r      = vars0.r;
    opts0.nIters = vars0.nIters;
    opts0.L0     = reshape(Xfft,[],nt);
    opts0.S0     = zeros(ny * nx,nt);
    [Lrpca, Srpca] = robustPCA(Y,vars0.lambdaL,vars0.lambdaS,opts0);
    
    % Construct LASSI initialization
    L0 = zeros(ny,nx,nt);
    S0 = Lrpca + Srpca;
elseif isfield(vars0,'mu1')
    % Add dependencies to path
    addpath('./deps_ktslr');
    
    % Run k-t SLR
    opts0.nItersO  = vars0.nItersO;
    opts0.nItersI  = vars0.nItersI;
    opts0.nItersCG = vars0.nItersCG;
    Xktslr = ktSLR(Y,A,vars0.p,vars0.mu1,vars0.mu2,Xfft,opts0);
    
    % Construct LASSI initialization
    L0 = zeros(ny,nx,nt);
    S0 = Xktslr;
else
    % Unrecognized algorithm
    error('Unrecognized algorithm');
end

% Run LASSI
opts.A       = A;
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
params.nLines = nLines;
params.seed   = seed; %#ok

% Save results
save(out,'Lhat','Shat','params','stats');
fprintf('DONE\n');

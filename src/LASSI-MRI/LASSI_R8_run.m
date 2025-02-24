function LASSI_R8_run(parFcn,idx)
% Syntax: LASSI_R8_run();
%         LASSI_R8_run(parFcn,idx);

% Get parameters
[vars, opts] = parFcn();

% Get parameters for this iteration
nR = numel(vars.r);
nL = numel(vars.lambdaL);
nS = numel(vars.lambdaS);
nB = numel(vars.lambdaB);
nD = numel(vars.dr);
nT = nR * nL * nS * nB * nD;
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
[ii, jj, kk, ll, mm] = ind2sub([nR, nL, nS, nB, nD],idx);
r        = vars.r(ii);
lambdaL  = vars.lambdaL(jj);
lambdaS  = vars.lambdaS(kk);
lambdaB  = vars.lambdaB(ll);
lambdaB0 = vars.lambdaB0(ll);
dr       = vars.dr(mm);

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

% Load undersampled data (Y, mask, samp, Xfft, Xtrue)
data = load(vars.inpath);
[ny, nx, nt] = size(data.Xtrue);
A = Emat_xyt(data.mask,data.samp,[ny, nx, nt]);

% Initialization
if isempty(vars.initpath)
    % FFT initialization
    init.Lhat = zeros(ny,nx,nt);
    init.Shat = data.Xfft;
else
    % Load initialization
    init = load(vars.initpath);
end

% Run LASSI
opts.A       = A;
if isfield(vars,'p')
    opts.p   = vars.p;
end
opts.r       = r;
opts.sdim    = [ny, nx, nt];
opts.dr      = dr;
opts.nIters  = vars.nIters;
opts.L0      = reshape(init.Lhat,[],nt);
opts.S0      = reshape(init.Shat,[],nt);
opts.Xtrue   = reshape(data.Xtrue,[],nt);
lB = logspace(log10(lambdaB0),log10(lambdaB),vars.nIters);
[Lhat, Shat, ~, ~, stats] = drpca(data.Y,lambdaL,lambdaS,lB,opts); %#ok

% Aggregate params
params.r        = r;
params.lambdaL  = lambdaL;
params.lambdaS  = lambdaS;
params.lambdaB  = lambdaB;
params.lambdaB0 = lambdaB0;
params.dr       = dr; %#ok

% Save results
save(out,'Lhat','Shat','params','stats');
fprintf('DONE\n');

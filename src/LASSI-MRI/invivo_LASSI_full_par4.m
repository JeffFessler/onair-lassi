function [vars, opts, vars0, opts0] = invivo_LASSI_full_par4()
% Syntax:   [vars, opts, vars0, opts0] = invivo_LASSI_full_par4();

% Knobs
vars.inpath   = 'invivo_full.mat';
vars.rawpath  = 'invivo_LASSI_full_par4/data.mat';
vars.outpath  = 'invivo_LASSI_full_par4.mat';
vars.nLines   = [5, 10, 15, 20, 25, 30];
vars.seed     = 1:5;
vars.r        = 4;
vars.lambdaL  = nan;
vars.lambdaS  = 0.005;
vars.lambdaB  = 0.003;
vars.lambdaB0 = 0.003;
vars.dr       = 5;
vars.nIters   = 50;

% DRPCA opts
%opts.A       = A;
%opts.sdim    = [ny, nx, nt];
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 2];
opts.type     = 'soft';
%opts.dr      = dr;
opts.ddim     = [prod(opts.pdim(1:(end - 1))), opts.pdim(end)];
%opts.nIters  = vars.nIters;
opts.nItersDB = 1;
opts.nItersLS = 5;
%opts.L0      = reshape(init.Lhat,[],nt);
%opts.S0      = reshape(init.Shat,[],nt);
opts.D0       = dctmtx(prod(opts.pdim));
%opts.Xtrue   = reshape(data.Xtrue,[],nt);
opts.accel    = false;
opts.tau      = 1;
opts.flag     = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Knobs
vars0.r       = nan;
vars0.lambdaL = 0.1;
vars0.lambdaS = 0.001;
vars0.nIters  = 250;

% RPCA opts
%opts0.A      = A;
%opts0.T      = T;
%opts0.r      = vars0.r;
%opts0.nIters = vars0.nIters;
%opts0.L0     = reshape(data.Xfft,[],nt);
%opts0.S0     = zeros(ny * nx,nt);
%opts0.Xtrue  = reshape(data.Xtrue,[],nt);
opts0.accel   = false;
opts0.tau     = 1;
opts0.flag    = 1;

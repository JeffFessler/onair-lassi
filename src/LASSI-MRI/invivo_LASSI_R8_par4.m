function [vars, opts] = invivo_LASSI_R8_par4()
% Syntax:   [vars, opts] = invivo_LASSI_R8_par4();

% Knobs
vars.inpath   = 'invivo_R8.mat';
vars.initpath = 'invivo_R8_lps_nrmse.mat';
vars.rawpath  = 'invivo_LASSI_R8_par4/data.mat';
vars.outpath  = 'invivo_LASSI_R8_par4.mat';
vars.r        = [1, 2, 3, 4, 5];
vars.lambdaL  = nan;
vars.lambdaS  = [0.005, 0.01, 0.02];
vars.lambdaB  = [0.0025, 0.005, 0.01];
vars.lambdaB0 = [0.0025, 0.005, 0.01];
vars.dr       = [1, 5];
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

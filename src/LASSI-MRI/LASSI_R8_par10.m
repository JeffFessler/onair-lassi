function [vars, opts] = LASSI_R8_par10()
% Syntax:   [vars, opts] = LASSI_R8_par10();

% Knobs
vars.inpath   = 'otazo_R8.mat';
vars.initpath = 'otazo_R8_lps_mse_0X.mat';
vars.rawpath  = 'LASSI_R8_par10/data.mat';
vars.outpath  = 'LASSI_R8_par10.mat';
vars.p        = 1; % SVT
vars.r        = nan;
vars.lambdaL  = linspace(0,2,11);
vars.lambdaS  = 0.02;
vars.lambdaB  = 0.03;
vars.lambdaB0 = 0.03;
vars.dr       = 1;
vars.nIters   = 50;

% DRPCA opts
%opts.A       = A;
%opts.sdim    = [ny, nx, nt];
opts.pdim     = [8, 8, 5];
opts.pgap     = [2, 2, 2];
opts.type     = 'hard';
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

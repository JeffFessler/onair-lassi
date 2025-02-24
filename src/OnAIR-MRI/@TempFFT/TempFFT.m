function T = TempFFT(dim,siz)
% Syntax:   T = TempFFT(dim);
%           T = TempFFT(dim,siz);

% Parse inputs
if ~exist('siz','var') || isempty(siz)
    % Default size
    siz = [];
end

% Create transform
T.adjoint = 0;
T.dim     = dim;
T.siz     = siz;
T         = class(T,'TempFFT');

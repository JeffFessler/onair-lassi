function  A = Emat_xyt(mask,b1,siz)
% Syntax:   A = Emat_xyt(mask,b1);
%           A = Emat_xyt(mask,b1,siz);

% Parse inputs
if ~exist('siz','var') || isempty(siz)
    % Default size
    siz = [];
end

% Create transform
A.adjoint = false;
A.mask    = mask;
A.b1      = permute(b1,[1 2 4 3]);
A.siz     = siz;
A         = class(A,'Emat_xyt');

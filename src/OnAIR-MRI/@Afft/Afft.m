function  A = Afft(mask,siz)
% Syntax:   A = Afft(mask);
%           A = Afft(mask,siz);

% Parse inputs
if ~exist('siz','var') || isempty(siz)
    % Default size
    siz = [];
end

% Create transform
A.adjoint = false;
A.mask    = mask;
A.siz     = siz;
A         = class(A,'Afft');

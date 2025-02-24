function D = formatDict(D)
% Syntax:   D = formatDict(D);

D = D - min(D(:));
D = D / max(D(:));

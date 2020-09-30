function Xr = tsvd(X, r)
% Syntax:   Xr = tsvd(X, r);

% Truncated SVD
[U, S, V] = svd(X, 'econ');
Xr = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)';

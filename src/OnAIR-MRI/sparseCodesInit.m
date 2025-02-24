function B0 = sparseCodesInit(D,X,pdim,pgap,sparsity)
% Syntax: B0 = sparseCodesInit(D,X,pdim,pgap,sparsity);

% Transformed patches
B0 = D' * X(patchInds(size(X),pdim,pgap));

% Threshold to desired sparisty
vals = sort(abs(B0(:)),'descend');
mu = vals(max(1,round((sparsity / 100) * numel(B0))));
B0 = hard(B0,mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hard thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = hard(X,mu)
X(abs(X) < mu) = 0;

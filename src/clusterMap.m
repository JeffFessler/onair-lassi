function Z = clusterMap(dim,y,Mp,Np)
%
% Syntax:       Z = clusterMap(dim,y,Mp,Np);
%               
% Inputs:       dim = size(A), where A is the underlying matrix from which
%               patches were extracted
%               
%               y is a 1 x size(Mp,2) vector of cluster assignments for
%               each patch of A
%               
%               [Mp, ~, Np] = patchInds(dim,...);
%               
% Outputs:      Z is a matrix of size(Z) = dim taking values in unique(y)
%               specifying the cluster membership of each pixel of A
%               
% Description:  Generates a majority vote cluster index map
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         April 10, 2016
%

% Parse inputs
n = size(Mp,2); % # of patches
N = max(Np(:)); % max # patches per pixels
d = prod(dim);  % # picels

% Cluster assignment matrix
Y = nan(d,N);
K = ones(d,1);
for j = 1:n
    Mpj    = Mp(:,j);
    KMpj   = K(Mpj);
    idx    = Mpj + d * (KMpj - 1);    
    Y(idx) = y(j);
    K(Mpj) = KMpj + 1;
end

% Majority vote
Z = reshape(mode(Y,2),dim);

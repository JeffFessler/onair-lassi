function [Mp, Dp, Np, Ip] = patchInds(adim,pdim,pgap,mask)
%
% Syntax:       Mp = patchInds(adim,pdim);
%               Mp = patchInds(adim,pdim,pgap);
%               Mp = patchInds(adim,pdim,pgap,mask);
%               [Mp, Dp, Np, Ip] = patchInds(adim,pdim);
%               [Mp, Dp, Np, Ip] = patchInds(adim,pdim,pgap);
%               [Mp, Dp, Np, Ip] = patchInds(adim,pdim,pgap,mask);
%               
% Inputs:       adim is a 1 x d vector specifying the dimensions of the
%               d-dimensional array, say A, you want to extract patches
%               from
%               
%               pdim is a 1 x d vector specifying the dimensions of the 
%               patches to extract
%               
%               [OPTIONAL] pgap is a 1 x d vector specifying the patch
%               strides (shifts) in each dimension. The default value is
%               pgap = pdim (non-overlapping patches)
%               
%               [OPTIONAL] mask is a logial array of size(M) = adim. If
%               specified, any patches containing indices (i1,...,id) for
%               which mask(i1,...,id) == false will be omitted
%               
% Outputs:      Mp is a p x n array whose columns specify the linear
%               indices of each patch. Here, p = prod(pdim), and n is the
%               total number of patches
%               
%               Dp is a 1 x d vector counting the number of patches along
%               each dimension of A. Note that if a mask is supplied and
%               patches are omitted, then n < prod(Dp)
%               
%               Np is an array of same size as A such that Np(i,j,...) is
%               the number of patches containing A(i,j,...)
%               
%               Ip is a 1 x n array containing the linear indices of the
%               patches (columns) of Mp in the Dp-dimensional grid from
%               which the patches were drawn. For example, if no mask is
%               supplied, Ip = 1:prod(Dp), while if a mask is specified,
%               the indcies corresponding to the omitted patches are absent
%               
% Description:  Generates an index matrix for fast multidimensional patch
%               extraction
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         April 10, 2016
%               December 3, 2016
%

% Parse inputs
if ~exist('pgap','var') || isempty(pgap)
    % Non-overlapping patches
    pgap = pdim;
end
if ~exist('mask','var')
    % No mask
    mask = [];
end
COMPUTE_NP = (nargout >= 3);

% Add trailing ones, if necessary
d    = max([numel(adim), numel(pdim), numel(pgap)]);
adim = [adim(:)', ones(1,d - numel(adim))];
pdim = [pdim(:)', ones(1,d - numel(pdim))];
pgap = [pgap(:)', ones(1,d - numel(pgap))];

% Compute base patch indices
I    = cell(1,d);
I{1} = 1:pdim(1);
for i = 2:d
    I{i} = 0:(pdim(i) - 1);
end
B = indmesh(adim,I{:});

% Compute offset indices
lidx = adim - pdim;
J    = cell(1,d);
for j = 1:d
    J{j} = [0:pgap(j):(lidx(j) - 1), lidx(j)];
end
O = indmesh(adim,J{:});

% Compute index matrix
Mp = bsxfun(@plus,B(:),O(:)');
if ~isempty(mask)
    % Apply mask
    Ip = find(~any(ismember(Mp,find(~mask)),1));
    Mp = Mp(:,Ip);
else
    Ip = 1:size(Mp,2);
end

% Compute patch array dimensions
Dp = cellfun(@numel,J);

% Compute count matrix, if requested
if COMPUTE_NP
    Np = zeros(adim);
    for i = 1:size(Mp,2)
        Np(Mp(:,i)) = Np(Mp(:,i)) + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate linear index mesh
function M = indmesh(dim,varargin)
n = numel(varargin);
M = varargin{1}(:);
for i = 2:n
    Mi = varargin{i}(:) * prod(dim(1:(i - 1)));
    M  = bsxfun(@plus,M,permute(Mi,circshift(1:i,[0 -1])));
end

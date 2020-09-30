function T = tilePatches(P,pdim,gdim,gap,val,scale)
%
% Syntax:       T = tilePatches(P,pdim,gdim);
%               T = tilePatches(P,pdim,gdim,gap);
%               T = tilePatches(P,pdim,gdim,gap,val);
%               T = tilePatches(P,pdim,gdim,gap,val,scale);
%               
% Inputs:       P is a p x n matrix whose columns contain image patches
%               
%               pdim is a 1 x 2 vector specifying the patch dimensions.
%               Note that p = prod(pdim)
%               
%               gdim is a 1 x 2 vector specifying the dimensions of the
%               patch grid to construct. Note that n = prod(gdim)
%               
%               [OPTIONAL] gap specifies the gap between patches in the
%               output tile image. The default value is gap = 1
%               
%               [OPTIONAL] val specifies the pixel value to fill the space
%               between tiles with. The default value is val = max(P(:))
%               
%               [OPTIONAL] scale = {true, false} determines whether to
%               scale each patch to the range [0, 1] before tiling. The
%               default value is scale = true
%               
% Outputs:      T is a matrix containing the tiled image patches from P
%               
% Description:  Generates a tiled image of the input patches
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         April 13, 2016
%               December 3, 2016
%

% Parse inputs
if ~exist('gap','var') || isempty(gap)
    % Default gap
    gap = 1;
end
if ~exist('val','var') || isempty(val)
    % Default val
    val = max(P(:));
end
if ~exist('scale','var') || isempty(scale)
    % Default scale
    scale = true;
end

% Scale columns
if scale
    % Scale P(:,j) to [0, 1]
    P = bsxfun(@minus,  P,min(P,[],1));
    P = bsxfun(@rdivide,P,max(P,[],1));
    P(isnan(P)) = 0;
end

% Tile columns
C = cell(gdim);
for j = 1:size(P,2)
    C{j} = reshape(P(:,j),pdim);
end
T = cell2tiles(C,gap,val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert cell to tile matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = cell2tiles(C,gap,val)
[m n] = size(C);
w     = sum(cellfun(@(c) size(c,2),C(1,:))) + (n - 1) * gap;
T     = catRow(C(1,:),gap,val);
HORZ  = val * ones(gap,w);
for i = 2:m
    T = cat(1,T,HORZ,catRow(C(i,:),gap,val));
end

% Concatenate tile row
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = catRow(R,gap,val)
T     = R{1};
VERT  = val * ones(size(T,1),gap);
for j = 2:numel(R)
    T = cat(2,T,VERT,R{j});
end

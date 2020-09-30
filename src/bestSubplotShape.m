function [nh, nw] = bestSubplotShape(n)
%
% Syntax:       [nh, nw] = bestSubplotShape(n);
%               
% Inputs:       n is the desired number of subplots
%               
% Outputs:      [nh, nw] are the recommended subplot dimensions
%               
% Description:  Recommends subplot dimensions for the given number of
%               subplots
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         February 2, 2017
%

% Get monitor aspect ratio
scrsz = get(0,'ScreenSize');
aspectRatio = scrsz(3) / scrsz(4);

% Compute optimal size
x = (1:n)';
A = bsxfun(@min,1 ./ x,aspectRatio ./ x');
A(x * x' < n) = 0;
[~, idx] = max(A(:));
[nh, nw] = ind2sub([n n],idx);

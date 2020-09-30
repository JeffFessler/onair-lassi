function fig = cfigure(dim,varargin)
%
% Syntax:       cfigure();
%               cfigure(dim);
%               cfigure(dim,...);
%               fig = cfigure();
%               fig = cfigure(dim);
%               fig = cfigure(dim,...);
%               
% Inputs:       [OPTIONAL] dim = [w, h] specifies the desired figure height
%               and width, in pixels. By default, the dimensions of
%               figure(...) are used
%               
%               [OPTIONAL] ... are valid arguments for figure(...)
%               
% Outputs:      fig is the figure handle
%               
% Description:  Resizes and centers the specified figure
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         April 12, 2016
%

% Get screen center
spos = get(0,'ScreenSize');
sxy0 = spos(1:2) + 0.5 * spos(3:4);

% Center figure
fig = figure(varargin{:});
if ~exist('dim','var') || isempty(dim)
    pos = get(fig,'Position');
    dim = pos(3:4);
end
set(fig,'Position',[(sxy0 - 0.5 * dim), dim]);

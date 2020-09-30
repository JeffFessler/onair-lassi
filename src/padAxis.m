function padAxis(ax,gap,xy)
%
% Syntax:       padAxis();
%               padAxis(ax);
%               padAxis(ax,gap);
%               padAxis(ax,gap,xy);
%               
% Inputs:       [OPTIONAL] ax is an axis handle or array of axis handles.
%               The default value is ax = gca()
%               
%               [OPTIONAL] gap > 0 is the desired relative x and y axis
%               paddings. The default value is gap = 0.05
%               
%               [OPTIONAL] xy = {'x','y','xy'} specifies which axes to pad.
%               The default value is xy = 'xy' (both axes)
%               
% Description:  Pads the current axis limits
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         June 28, 2016
%

% Parse inputs
if ~exist('ax','var') || isempty(ax)
    % Default axis
    ax = gca;
end
if ~exist('gap','var') || isempty(gap)
    % Default gap
    gap = 0.05;
end
if ~exist('xy','var') || isempty(xy)
    % Default xy
    xy = 'xy';
end

% Pad axis limits
for i = 1:numel(ax)
    for j = 1:numel(xy)
        padAxisLimits(ax(i),gap,xy(j));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pad axis limits
function padAxisLimits(ax,gap,xy)
isLog  = strcmpi(get(ax,[xy, 'scale']),'log');
oldLim = get(ax,[xy, 'lim']);
if isLog
    oldLim = log(oldLim);
end
newLim = oldLim + gap * diff(oldLim) * [-1, 1];
if isLog
    newLim = exp(newLim);
end
set(ax,[xy, 'lim'],newLim);

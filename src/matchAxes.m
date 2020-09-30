function matchAxes(arg1,gap,xy)
%
% Syntax:       matchAxes();
%               matchAxes(fig);
%               matchAxes(ax);
%               matchAxes(fig,gap);
%               matchAxes(ax,gap);
%               matchAxes(fig,gap,xy);
%               matchAxes(ax,gap,xy);
%               
% Inputs:       [OPTIONAL] fig is a figure handle. By default, all axes of
%               gcf() are matched
%               
%               [OPTIONAL] ax is an axis handle or array of axis handles.
%               By default, all axes of gcf() are matched
%               
%               [OPTIONAL] gap > 0 is the desired relative axis padding.
%               The default value is gap = 0.05
%               
%               [OPTIONAL] xy = {'x','y','xy'} specifies which axes to
%               format. The default value is xy = 'xy' (both axes)
%               
% Descprition:  Match the axis limits of the given axes
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         June 28, 2016
%

% Parse inputs
if ~exist('arg1','var') || isempty(arg1)
    arg1 = gcf();
end
if ~exist('gap','var') || isempty(gap)
    % Default padding
    gap = 0.05;
end
if ~exist('xy','var') || isempty(xy)
    % Default axes
    xy = 'xy';
end

% Get axes
if strcmpi(get(arg1(1),'type'),'figure')
    % A figure handle was provided
    fig = arg1;
    ax  = setdiff(findobj(fig,'type','axes'), ...
                  findobj(fig,'tag','legend'));
else
    % Axis handle(s) were provided
    ax = arg1;
end

% Set tight limits
for i = 1:numel(xy)
    setTightLimits(ax,xy(i));    
end

% Pad axes
padAxis(ax,gap,xy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set tight limits for given axes
function setTightLimits(ax,xy)
% Get axis limits
n = numel(ax);
lim = nan(n,2);
for i = 1:n
    lim(i,:) = get(ax(i),[xy, 'lim']);
end

% Set tight limits
tlim = [min(lim(:,1)), max(lim(:,2))];
set(ax,[xy, 'lim'],tlim);

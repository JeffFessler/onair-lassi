function varargout = deps(path,kind,verbose,type,varargin)
%
% Syntax:       deps();
%               deps(path);
%               deps(path,kind);
%               deps(path,kind,verbose);
%               deps(path,kind,verbose,type);
%               [files, M] = deps(...);
%               
% Inputs:       [OPTIONAL] path can be one of the following:
%               
%                   (a) The relative, absolute, or partial path to an .m
%                       file or a directory
%                       
%                   (b) A string like '*.m' such that dir(path) returns a
%                       list of .m files. Directories and other non m-files
%                       returned by dir(path) are ignored
%                   
%               The default value is path = '.' (current directory)
%               
%               [OPTIONAL] kind specifies what type of function
%               dependencies to list. The options are:
%               
%                   kind = {'user'   , User-defined functions
%                          {'builtin', MATLAB built-in/toolbox functions
%                          {'all'    , All (non-Java) functions
%               
%               The default value is kind = 'user'
%               
%               [OPTIONAL] verbose = [0 1] determines what results to
%               print. The options are:
%               
%                   0: No printing
%                   1: Print dependencies (excluding input files)
%               
%               The default value is verbose = {0, if nargout >= 1
%                                              {1, if nargout == 0
%               
%               [OPTIONAL] type = [0 1 2 3] determines what type of
%               dependency matrix to display. The options are:
%               
%                   0: No display
%                   1: Full dependency matrix
%                   2: Same as 1 with orphan callers/called removed
%                   3: Same as 2 with input files removed
%               
%               The default value is type = 0
%               
% Outputs:      files is an nFiles x 1 cell array containing the absolute
%               path to each dependent file (including input files)
%               
%               M is an nFiles x nFiles sparse dependency matrix such that
%               
%                   M(i,j) = {1, if files{i} calls files{j}
%                            {0, otherwise
%               
% Description:  Lists all dependencies of the given .m file or each .m file
%               in the given directory
%               
% Examples:     % List dependencies of this file
%               deps('deps.m');
%               
%               % Display dependency matrix for current directory
%               deps([],[],[],1);
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         June 28, 2016
%               August 6, 2016
%

% Handle [?] callback
if ~isempty(varargin)
    % Print callers
    printCallers(varargin{:});
    return;
end

% Parse inputs
if ~exist('path','var') || isempty(path)
    % Default path
    path = '.';
end
if ~exist('kind','var') || isempty(kind)
    % Default kind
    kind = 'user';
end
if ~exist('verbose','var') || isempty(verbose)
    % Default verbosity
    verbose = double(~nargout);
end
if ~exist('type','var') || isempty(type)
    % Default type
    type = 0;
end

% Initial file list
switch exist(path,'file')
    case 2
        % File
        files = {stdPath(which(path))};
    case 7
        % Directory
        W     = what(path);
        names = W.m;
        files = absolutePath(getBaseDir(path),names(:));
    otherwise
        % Search pattern
        D     = dir(path);
        filt  = @(d) ~d.isdir && any(regexp(d.name,'.m$'));
        names = {D(arrayfun(filt,D)).name};
        files = absolutePath(getBaseDir(path),names(:));
        
        % Invalid path
        if isempty(files)
            error('No files/directories matching "%s"',path);
        end
end
n0 = numel(files);

% Get dependencies
idx = 0;
Mi  = [];
Mj  = [];
while idx < numel(files)
    idx   = idx + 1;
    calls = getFunctionCalls(files{idx},kind);
    files = [files; calls(~ismember(calls,files))]; %#ok
    Mi    = [Mi; repmat(idx,numel(calls),1)];  %#ok
    Mj    = [Mj; find(ismember(files,calls))]; %#ok
end

% Construct dependency matrix
nFiles = numel(files);
M      = sparse(Mi,Mj,1,nFiles,nFiles);

% Print dependencies, if requested
if verbose > 0
    % Setup caller checking
    assignin('base','deps_FILES',files);
    assignin('base','deps_M',M);
    
    % Print files
    printOrganizedFiles(files((n0 + 1):nFiles));
end

% Display dependency matrix, if requested
if type == 1
    % Full matrix
    displayDependencyMatrix(path,M,files,files,n0,n0);
elseif type > 1
    % Compact dependency matrix
    c = any(M,1);
    if type == 3
        % Remove input files from callee list
        c(1:n0) = false;
    end
    r = any(M(:,c),2);
    filesy = files(r);
    filesx = files(c);
    y0 = nnz(r(1:n0));
    x0 = nnz(c(1:n0));
    displayDependencyMatrix(path,M(r,c),filesy,filesx,y0,x0);
end

% Populate outputs
varargout = populateOutputs(nargout,files,M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Print dependencies that call given file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printCallers(file,files,M)
callers = files(any(M(:,ismember(files,file)),2));
if ~isempty(callers)
    [~, name, ext] = fileparts(file);
    fprintf('<strong>Functions that call "%s%s"</strong>\n',name,ext);
    printOrganizedFiles(callers);
end

% Construct absolute path(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function apath = absolutePath(base,names)
absFcn = @(n) sprintf('%s/%s',base,n);
if ~iscell(names)
    apath = absFcn(names);
else
    apath = cellfun(absFcn,names,'UniformOutput',false);
end

% Get base directory for the given path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function base = getBaseDir(path)
base = fullpath(path);                  % Handles files on MATLAB path
if isempty(base)
    base = fullpath(fileparts(path));   % Handles relative search patterns
    if isempty(base)
        base = stdPath(pwd());          % Handles search patterns in '.'
    end
end

% Get paths to functions of the given kind called by the given function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calls = getFunctionCalls(base,kind)
info = getcallinfo(base,'flat');
fcns = unique([info.calls.fcnCalls.names(:);    % Standard function calls
               info.calls.dotCalls.names(:);    % "Dot" function calls
               info.calls.atCalls.names(:)]);   % "At" function calls
if isempty(fcns)
    % No function calls found
    calls = {};
    return;
end
mPath = @(f) getMfilePath(f,base,kind);
calls = cellfun(mPath,fcns,'UniformOutput',false);
calls = unique(calls(~cellfun(@isempty,calls)));

% Return path to .m file of given kind, or empty if it does not exist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function path = getMfilePath(func,base,kind)
path   = stdPath(which([func, '.m'],'in',base));
mlroot = stdPath(matlabroot());
if any(strfind(path,'is a Java method'))
    % Java method
    path = '';
elseif strncmp(path,'built-in',8)
    % Built-in method
    if strncmpi(kind,'user',1)
        path = '';
    else
        % Extract path
        path = regexprep(path,'.*\((.+)\)','$1');
    end
elseif strncmp(path,mlroot,numel(mlroot))
    % Toolbox (or built-in?) method)
    if strncmpi(kind,'user',1)
        path = '';
    end
elseif strncmpi(kind,'builtin',1)
    % User-defined function
    path = '';
end

% Print organized files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printOrganizedFiles(files)
paths = cellfun(@fileparts,files,'UniformOutput',false);
[paths, ~, inds] = unique(paths);   % Organized by directory
rootdir = root();                   % MATLAB directory
nPaths  = numel(paths);
for i = 1:nPaths
    printPath(paths{i},rootdir);    % Print directory
    idx      = find(inds == i);
    nMatches = numel(idx);
    for j = 1:nMatches
        printMatch(files{idx(j)});  % Print file
    end
end
if nPaths > 0
    fprintf('\n');
end

% Print path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printPath(fpath,rootdir)
cmd   = sprintf('go(''%s'');',fpath);
link  = sprintf('<a href="matlab:%s">[GO]</a>',cmd);
fpath = regexprep(fpath,['^', rootdir],'~');
fprintf('\n%s: %s\n',link,fpath);

% Print match
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printMatch(fpath)
[~, name, ext] = fileparts(fpath);
cmde  = sprintf('edit(''%s'');',fpath);
linke = sprintf('<a href="matlab:%s"><-</a>',cmde);
cmdc  = sprintf('deps([],[],[],[],''%s'',deps_FILES,deps_M);',fpath);
linkc = sprintf('<a href="matlab:%s">[?]</a>',cmdc);
str   = sprintf('%s%s',name,ext);
fprintf(' %s : %s %s\n',linke,str,linkc);

% Get base filename
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function name = getBaseFilename(file)
[~, name] = fileparts(file);

% Standardize path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function path = stdPath(path)
%path = regexprep(path,'\\|/',filesep());   % Use platform's file separator
path  = regexprep(path,'\\','/');           % Use Unix file separator

% Populate outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function argout = populateOutputs(narg,varargin)
argout = varargin(1:narg);

% Display dependency matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayDependencyMatrix(path,M,filesy,filesx,y0,x0)
% Constants
THETA  = 90;                        % x-label rotation, in degrees
ALPHA  = 0.15;                      % Patch transparency
COLORS = [000, 147, 059;            % Google green
          249, 001, 001;            % Google red
          030, 030, 030;            % Soft black
          002, 102, 200] / 255;     % Google blue
FONT   = {'FontUnits','pixels', ...
          'Interpreter','None'};    % Font options

% Parse inputs
[ny, nx] = size(M);
fmtName  = @(n) sprintf(' %s ',getBaseFilename(n));
namesy   = cellfun(fmtName,filesy,'UniformOutput',false);
namesx   = cellfun(fmtName,filesx,'UniformOutput',false);
[i, j]   = find(M);
idx      = [(i <= y0) & (j <= x0), ...
            (i <= y0) & (j >  x0), ...
            (i >  y0) & (j <= x0), ...
            (i >  y0) & (j >  x0)];

% Display dependency graph
fig = figure('name','deps', ...
             'MenuBar','none', ...
             'DockControl','off', ...
             'NumberTitle','off');
ax(1) = axes('NextPlot','add');
for k = 1:4
    plot(j(idx(:,k)),i(idx(:,k)),'.', ...
                                 'Parent',ax(1), ...
                                 'MarkerSize',20, ...
                                 'Color',COLORS(k,:));
end
set(ax(1),'XLim',[0, nx + 1], ...
          'YLim',[0, ny + 1], ...
          'YDir','reverse', ...
          'Visible','off');

% Display background image
y   = 2 * ((1:ny)' > y0);
x   = 1 * ((1:nx)' > x0);
P   = 1 + bsxfun(@plus,y,x');
ax(2) = axes('Position',get(ax(1),'Position'));
image(P,'Parent',ax(2));
colormap((1 - ALPHA) + ALPHA * COLORS);
box(ax(2),'on');
grid(ax(2),'on');
set(ax(2),'XLim',[0, nx + 1], ...
          'YLim',[0, ny + 1], ...
          'YDir','reverse', ...
          'XTick',[], ...
          'YTick',[]);
titleStr = sprintf('Dependency graph for "%s"',path);
th = title(titleStr,'Parent',ax(2),FONT{:});
fontSize = get(th,'FontSize');

% Format figure
uistack(ax(2),'bottom');
uistack(ax(1),'top');
set(fig,'Units','pixels');
set(ax,'Units','pixels');

% Add y-labels
xOff = min(get(ax(2),'XLim'));
set(ax(2),'YTick',0:ny);
set(ax(2),'YTickLabel',[]);
ty = text(repmat(xOff,1,ny),1:ny,namesy,'HorizontalAlignment','right', ...
                                        FONT{:});
ty(ny + 1) = text(xOff,0,' CALLER ','HorizontalAlignment','right', ...
                                    'FontWeight','Bold', ...
                                    FONT{:});

% Add x-labels
yOff = max(get(ax(2),'YLim'));
set(ax(2),'XTick',0:nx);
set(ax(2),'XTickLabel',[]);
tx = text(1:nx,repmat(yOff,1,nx),namesx,'HorizontalAlignment','right', ...
                                        'Rotation',THETA, ...
                                        FONT{:});
tx(nx + 1) = text(0,yOff,' CALLED ','HorizontalAlignment','right', ...
                                    'Rotation',THETA, ...
                                    'FontWeight','Bold', ...
                                    FONT{:});

% Resize figure
kappa  = (ny + 2) / (nx + 2);
gap0   = [getMaxDim(ty) + 5, getMaxDim(tx)];
gap1   = 2.5 * fontSize * [1, 1];
gap    = gap0 + gap1;
height = (ny + 1) * fontSize + gap(1);
resizeFcn = @(s,e) resizeAxis(fig,ax,gap0,gap1);
set(fig,'ResizeFcn',resizeFcn);
resizeFcn();
setAspectRatio(fig,kappa,gap,height);

% Add menus to figure
filem = uimenu(fig,'Label','File');
uimenu(filem,'Label','Save...', ...
             'Callback',@(s,e) saveGraph(fig), ...
             'Accelerator','S');
uimenu(filem,'Label','Close', ...
             'Callback',@(s,e) delete(fig), ...
             'Accelerator','W');
viewm = uimenu(fig,'Label','View');
uimenu(viewm,'Label','Match Aspect Ratio', ...
             'Callback',@(s,e) setAspectRatio(fig,kappa,gap), ...
             'Accelerator','A');

% Get max text box dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = getMaxDim(t)
set(t,'Units','pixels');
w = max(cellfun(@(e) max(e(3:4)),get(t,'Extent')));
set(t,'Units','data');

% Save graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveGraph(fig)
% Get output path from user
outpath = inputdlg('Output file:', ...
                   'Save dependency graph',[1 50], ...
                   {'./deps.pdf'});
drawnow();
if isempty(outpath)
    % User didn't want to save after all...
    return;
end
outpath = outpath{1};

% Write file
[~, ~, ext] = fileparts(outpath);
export_fig(fig,['-', ext(2:end)],'-transparent','-nocrop',outpath);
fprintf('File ''%s'' written\n',outpath);

% Set figure aspect ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setAspectRatio(fig,kappa,gap,height)
pos = get(fig,'Position');
if nargin == 4
    % Set height
    pos(4) = height;
end
pos(3) = (pos(4) - gap(2)) / kappa + gap(1); % Aspect-perserving width
set(fig,'Position',pos);

% Resize axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resizeAxis(fig,ax,gap0,gap1)
fpos = get(fig,'Position');
set(ax,'Position',[gap0, fpos(3:4) - (gap0 + gap1)]);

function paths = findMatchingFiles(basepath)
%
% Syntax:       paths = findMatchingFiles(basepath);
%               
% Inputs:       basepath is the base filename (absolute or relative) to
%               search
%               
% Outputs:      paths is a cell array of paths to files matching the
%               regex [path, name, '*', ext], where [path, name, ext] =
%               fileparts(basepath), sorted by appended number (if any)
%               
% Description:  Return all files with the given base filename, sorted by
%               appended number (if any)
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         January 26, 2017
%

% Parse base path
basepath = regexprep(basepath,'\','/');
[path, name, ext] = fileparts(basepath);
if ~isempty(path)
    path = [path, '/'];
end

% Find matching files
list = dir([path, name, '*', ext]);
names = {list.name};
names = names(:);

% Sort by appended number (if any)
[~, idx] = sort(cellfun(@getNumber,names));
names = names(idx);

% Return fullpaths to matching files
fcn = @(name) [path, name];
paths = cellfun(fcn,names,'UniformOutput',false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = getNumber(name)
[~, str, ~] = fileparts(name);
n = str2double(regexp(str,'\d+$','match'));
if isempty(n)
    n = nan;
end

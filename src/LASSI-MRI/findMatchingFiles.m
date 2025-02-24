function paths = findMatchingFiles(basepath)
% Syntax:       paths = findMatchingFiles(basepath);
% Description:  Return all files with the given base filename, sorted by
%               appended number

% Parse base filepath
basepath = regexprep(basepath,'\','/');
[path, name, ext] = fileparts(basepath);
if ~isempty(path)
    path = [path, '/'];
end

% Find matching filenames
list = dir([path, name, '*', ext]);
names = {list.name};
names = names(:);

% Sort by appended number
[~, idx] = sort(cellfun(@getNumber,names));
names    = names(idx);

% Return complete paths to matching files
fcn   = @(name) [path, name];
paths = cellfun(fcn,names,'UniformOutput',false);

%
% Local functions
%

% Get number
function n = getNumber(name)
[~,str,~] = fileparts(name);
n = str2double(regexp(str,'\d+$','match'));
if isempty(n)
    n = nan;
end

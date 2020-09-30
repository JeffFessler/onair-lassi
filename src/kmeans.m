function [C, y] = kmeans(X,k,mode,maxIters,flag)
%
% Syntax:       C = kmeans(X,k);
%               C = kmeans(X,k,'++');
%               C = kmeans(X,k,'++',maxIters);
%               C = kmeans(X,k,'++',maxIters,flag);
%               [C, y] = kmeans(X,k);
%               [C, y] = kmeans(X,k,'++');
%               [C, y] = kmeans(X,k,'++',maxIters);
%               [C, y] = kmeans(X,k,'++',maxIters,flag);
%               
% Inputs:       X is a d x n matrix containing n d-dimensional samples
%               
%               k is the desired number of clusters
%               
%               [OPTIONAL] '++' chooses initial cluster centers using the
%               kmeans++ algorithm. By default, the initial clusters are
%               chosen randomly
%               
%               [OPTIONAL] maxIters is the maximum allowed number of
%               iterations. The default value is maxIters = 100
%               
%               [OPTIONAL] flag = {true, false} determines whether to print
%               iteration stats to the command window. The default value is
%               flag = false
%               
% Outputs:      C is a d x k matrix whose columns contain the centers of 
%               the k clusters
%               
%               y is a 1 x n vector of group labels y(i) = {1,...k} of each
%               input sample X(:,i)
%               
% Description:  Clusters the input data using the k-means algorithm
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         April 10, 2016
%

% Parse inputs
if ~exist('mode','var') || isempty(mode)
    % Default mode
    mode = '';
end
if ~exist('maxIters','var') || isempty(maxIters)
    % Default max # iterations
    maxIters = 100;
end
if ~exist('flag','var') || isempty(flag)
    % Default flag
    flag = false;
end

% Initialize printing
if flag
    % Stats-printing function
    iterFmt = sprintf('%%0%dd',ceil(log10(maxIters + 1)));
    out     = printFcn('Iter',iterFmt ,'time','%.2fs');
    if strcmpi(mode,'++')
        % kmeans++
        fprintf('***** kmeans++ *****\n');
    else
        % Standard kmeans
        fprintf('***** kmeans *****\n');
    end
end

% Initialize cluster centers
[d, n] = size(X);
if strcmpi(mode,'++')
    % k-means++ centers
    C      = zeros(d,k);
    idx    = randi(n);
    C(:,1) = X(:,idx);
    inds   = setdiff(1:n,idx);
    for i = 1:(k - 1)
        Di  = min(pdist2(C(:,1:i)',X(:,inds)').^2,[],1);
        Fi  = cumsum(Di) / sum(Di);
        idx = find(rand() <= Fi,1,'first');
        C(:,i + 1) = X(:,inds(idx));
        inds(idx)  = [];
    end
else
    % Random centers
    C = X(:,randperm(n,k));
end

% Perform clustering
y    = nan(1,n);
it   = 0;
done = false;
while ~done
    % Initialize iteration
    itimer = tic;
    it     = it + 1;
    yLast  = y;
    
    % Update clusters
    [~, y] = min(pdist2(X',C')',[],1); %#ok
    
    % Update centroids
    for i = 1:k
        C(:,i) = mean(X(:,y == i),2);
    end
    
    % Check convergence
    done = all(y == yLast) || (it >= maxIters);
    
    % Print status
    if flag
        out(it,toc(itimer));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Print function
function out = printFcn(varargin)
str = [sprintf('%s[%s] ',varargin{:}), '\n'];
out = @(varargin) fprintf(str,varargin{:});

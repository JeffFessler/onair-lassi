function X = interpVideo(X,M,algo,method)
%
% Syntax:       X = interpVideo(X,M);
%               X = interpVideo(X,M,algo,method);
%               
% Inputs:       X is an ny x nx x nt video
%               
%               M is an ny x nx x nt observation mask
%               
%               [OPTIONAL] {algo,method} specify the interpolation strategy
%               to use. They can take values in the table below. The
%               default values are {'paint2',1}
%               
%                    algo | method
%               ----------+--------------------
%                'paint2' | 0,1,2,3,4
%                'paint3' | 1
%                 'grid2' | 'cubic','linear'
%                  'tri2' | 'natural','linear'
%                  'tri3' | 'natural','linear'
%               
% Outputs:      X is a copy of the input video with missing values
%               interpolated
%               
% Description:  Interpolates the values of X(i,j,k) for which M(i,j,k) == 0
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         January 26, 2017
%

% Parse inputs
if ~exist('algo','var') || ~exist('method','var')
    % Default method
    algo = 'paint2';
    method = 1;
end

% Perform interpolation
switch lower(algo)
    case 'paint2'
        % 2D interpolation via inpaint_nans
        Y = X;
        Y(~M) = nan;
        for k = 1:size(X,3)
            X(:,:,k) = inpaint_nans(Y(:,:,k),method);
        end
    case 'paint3'
        % 3D interpolation via inpaint_nans3
        Y = X;
        Y(~M) = nan;
        X = inpaint_nans3(Y,method);
    case 'grid2'
        % 2D interpolation via griddata
        for k = 1:size(X,3)
            Xk = X(:,:,k);
            Mk = M(:,:,k);
            [i1, j1] = find(Mk);
            [i0, j0] = find(~Mk);
            Xk(~Mk) = griddata(i1,j1,Xk(Mk),i0,j0,method); %#ok
            X(:,:,k) = Xk;
        end
    case 'tri2'
        % 2D interpolation via TriScatteredInterp()
        for k = 1:size(X,3)
            Xk = X(:,:,k);
            Mk = M(:,:,k);
            [i1, j1] = find(Mk);
            [i0, j0] = find(~Mk);
            F = TriScatteredInterp(i1,j1,Xk(Mk),method);
            Xk(~Mk) = F(i0,j0);
            X(:,:,k) = Xk;
        end
    case 'tri3'
        % 3D interpolation via TriScatteredInterp()
        [i1, j1, k1] = ind2sub(size(M),find(M));
        [i0, j0, k0] = ind2sub(size(M),find(~M));
        F = TriScatteredInterp(i1,j1,k1,X(M),method);
        X(~M) = F(i0,j0,k0);
    otherwise
        % algo not supported
        error('algo ''%s'' not supported',algo);
end

% Fill nans
X(isnan(X)) = 0;

function Y = patchRecon(Z,Mp,Np)
%
% Syntax:       Y = patchRecon(Z,Mp,Np);
%               
% Inputs:       Z is a p x n matrix of patch data
%               
%               Mp and Np are the outputs of [Mp, ~, Np] = patchInds(...)
%               
% Outputs:      Y is an array of size(Np) containing the reconstructed
%               patch data
%               
% Description:  Reconstructs a dataset from its (possibly overlapping) 
%               patches. For example, 
%               
% Example:      One always has the relationship
%                   
%                   Y == patchRecon(Y(Mp),Np)
%               
%               where
%               
%                   [Mp, ~, Np] = patchInds(size(Y),...);
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         January 22, 2017
%               March 4, 2017
%

% Parse inputs
n = size(Mp,2); % Number of patches
dim = size(Np); % Output image size

% Patch reconstruction
Y = zeros(dim);
for i = 1:n
    Mpi = Mp(:,i);
    Y(Mpi) = Y(Mpi) + Z(:,i); % Pixel-wise accumulation
end
Y = Y ./ max(1,Np); % Pixel-wise average

function Y = corrupt(X,SNR,type)
%
% Syntax:       Y = corrupt(X,SNR);
%               Y = corrupt(X,SNR,type);
%               
% Inputs:       X is an arbitrary array
%               
%               SNR is the desired signal-to-noise ratio, in dB
%               
%               [OPTIONAL] type = {'gaussian','poisson'} specifies the type
%               of noise to add. The default value is type = 'gaussian'
%               
% Outputs:      Y is a copy of X corrupted by noise of the given type and
%               SNR
%               
% Description:  Corrupts data with noise of given type and SNR
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         December 15, 2016
%               January 23, 2017
%

% Parse inputs
if ~exist('type','var') || isempty(type)
    % Default noise type
    type = 'gaussian';
end

% Add noise
if strncmpi(type,'gaussian',1)
    % Gaussian noise
    if isreal(X)
        % Real-valued
        sigma = 10^(-SNR / 20) * norm(X(:)) / sqrt(numel(X));
        Y = X + sigma * randn(size(X));
    else
        % Complex-valued
        sigma = 10^(-SNR / 20) * norm(X(:)) / sqrt(2 * numel(X));
        Y = X + sigma * (randn(size(X)) + 1i * randn(size(X)));
    end
elseif strncmpi(type,'poisson',1)
    % Poisson noise
    alpha = (sum(X(:)) / sum(X(:).^2)) * 10^(0.1 * SNR);
    Y = poissrnd(alpha * X) / alpha;
else
    % Unsupported type
    error('Noise type "%s" is not supported',type);
end

%{
% Validate SNR
SNRact = 10 * log10(sum(X(:).^2) / sum((Y(:) - X(:)).^2)) %#ok
%}

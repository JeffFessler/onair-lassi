function [Y, Xtrue, X0, b1] = generateCardiacPerfData2(p,SNR,data)
% Syntax:   [Y, Xtrue, X0, b1] = generateCardiacPerfData2(p,SNR,data);

% Load ground truth data
if ischar(data)
    % Path supplied, so load data
    data = load(data);
end
Xtrue = data.Xtrue;
b1 = data.samp; % Coil sensitivities

% Undersample data
[Y, mask] = subsampleData(data.Ytrue,p,data.ps);

% Add noise, if requested
if ~isempty(SNR) && isfinite(SNR)
    inds      = repmat(mask,[1 1 1 data.nc]);
    np        = nnz(inds);
    SNR2sigma = @(SNR,Y) exp(-SNR / 20) * (norm(Y(:)) / sqrt(numel(Y))) / sqrt(2);
    Y(inds)   = Y(inds) + SNR2sigma(SNR,Y(inds)) * (randn(np,1) + 1i * randn(np,1));
end

%--------------------------------------------------------------------------
% Compute initial iterate
%--------------------------------------------------------------------------
%{
% Method #1: Back-projection
X0 = A' * Y;
%}

%{
% Method #2: Data-sharing back propagation
%
% BUGGY???
%
samp   = flipdim(flipdim(data.samp,1),2);
Yfill  = data_share_fill(reshape(Y(abs(Y) ~= 0),[],data.nc),mask);
center = @(X) fftshift(fftshift(X,1),2);
xc     = permute(center(ifft2(center(Yfill))),[1 2 4 3]);
Xfft   = squeeze(sum(xc .* conj(repmat(samp,[1 1 1 data.nt])),3)) ./ ...
         repmat(sum(abs(samp).^2,3),[1 1 data.nt]);
X0     = sqrt(data.nx * data.ny) * flipdim(flipdim(Xfft,1),2);
%}

% Method #2: Data-sharing back propagation
siz = size(mask);
Yfill = data_share_fill(reshape(Y(abs(Y) ~= 0),[],data.nc),mask);
Afill = Emat_xyt(1,b1,siz);
X0 = reshape(Afill' * Yfill,siz);

%{
% Method #3: Least-squares
X0 = A' * Y;
tau = 1.5 / normA^2; % Numerator in (0,2)
Niters = 20;
for i = 1:Niters
    X0last = X0;
    X0 = X0 - tau * (A' * (A * X0 - Y));
    delta = norm(X0(:) - X0last(:)) / norm(X0last(:));
    fprintf('Iteration %i/%i (delta = %.3g)\n',i,Niters,delta);
end
%}

% Scale initial iterate
alpha = fminsearch(@(K) norm(abs(Xtrue(:)) - abs(K * X0(:))),1);
X0 = alpha * X0;
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subsample data
function [Ys mask] = subsampleData(Y,p,ps)

% Parse inputs
[ny nx nt nc] = size(Y);
num = round(p * nx);

% Undersample data
Ys = zeros(ny,nx,nt,nc);
mask = false(ny,nx,nt);
for i = 1:nt
    cols = sampleCols(num,ps);
    mask(:,cols,i) = true;
    Ys(:,cols,i,:) = Y(:,cols,i,:);
end

% Sample columns
function inds = sampleCols(n,ps)

% Sample according to given distribution
inds = [];
while numel(inds) < n
    P       = cumsum(ps);
    [~,idx] = histc(rand(1,n),[0, P / P(end)]);
    idx     = unique(idx);
    ps(idx) = 0;
    inds    = [inds, idx(randperm(numel(idx)))]; %#ok
end
inds = sort(inds(1:n));

% Data-sharing (0th order) k-space interpolation
function Yfull = data_share_fill(Y,mask)
% Syntax:   Yfull = data_share_fill(Y,mask);

Nx = size(mask,1);
Ny = size(mask,2);
Nf = size(mask,3);
Ns = size(Y,1);
Nc = size(Y,2); 
flip_samp = flipdim(mask,3);
samps = {mask,flip_samp};
zfill = reshape(embed(Y,logical(mask)),Nx,Ny,Nf,Nc);
for ii = 1:2
	curr_samp = samps{ii};
	%interp_ind = cumsum(int32(curr_samp),3);
	interp_ind = cumsum((curr_samp),3);
	%interp_ind_pos = interp_ind + int32(interp_ind == 0);
	interp_ind_pos = interp_ind + (interp_ind == 0);
	siip(:,:,:,ii) = interp_ind_pos; % keep!!
    
	%full_interp_ind_pos = Nx*Ny*(interp_ind_pos(:)-1) + int32(repmat((0:Nx*Ny-1)', [Nf 1])) + 1;
	full_interp_ind_pos = Nx*Ny*(interp_ind_pos(:)-1) + (repmat((0:Nx*Ny-1)', [Nf 1])) + 1;
    
	if Nc > 1
		%full_interp_ind_pos = repmat(full_interp_ind_pos,[Nc 1]) + int32(kron((0:Nc-1)',Nx*Ny*Nf*ones(Nx*Ny*Nf,1)));	
		full_interp_ind_pos = repmat(full_interp_ind_pos,[Nc 1]) + (kron((0:Nc-1)',Nx*Ny*Nf*ones(Nx*Ny*Nf,1)));	
    end
    
	[~, inds] = sort(curr_samp,3,'descend');
	sinds(:,:,:,ii) = inds; % keep!!
    
	full_inds = Nx*Ny*(inds(:)-1) + repmat((1:Nx*Ny)',[Nf 1]);
	if Nc > 1
		full_inds = repmat(full_inds,[Nc 1]) + kron((0:Nc-1)',Nx*Ny*Nf*ones(Nx*Ny*Nf,1));
    end
    
	if ii == 2
		curr_zfill = vec(flipdim(zfill,3));
	else
		curr_zfill = vec(zfill);
	end
	nonzero_elms(:,ii) = curr_zfill(full_inds);
	zero_less_fill = reshape(nonzero_elms(:,ii),Nx,Ny,Nf,Nc);
	curr_full_data = reshape(zero_less_fill(full_interp_ind_pos),Nx,Ny,Nf,Nc);
	if ii == 2 
		curr_full_data  = flipdim(curr_full_data ,3);
	end
	sfull_data(:,:,:,:,ii) = curr_full_data;
end
siip1 = double(siip(:,:,:,1));
c1_ndx = (siip1-1)*Nx*Ny+repmat(reshape(0:Nx*Ny-1,Nx,Ny),[1 1 Nf])+1;
frames1 = vec(sinds(:,:,:,1));
c1 = vec(reshape(frames1(c1_ndx),Nx,Ny,Nf));
siip2 = double(flipdim(siip(:,:,:,2),3));
siip2 = repmat(siip2(:,:,1) + 1,[1 1 Nf]) - siip2;
c2_ndx = (siip2-1)*Nx*Ny+repmat(reshape(0:Nx*Ny-1,Nx,Ny),[1 1 Nf])+1;
frames2 = vec(sinds(:,:,:,1)); % need to have 1!
c2 = vec(reshape(frames2(c2_ndx),Nx,Ny,Nf));
c3 = vec(reshape(kron((1:Nf),ones(Nx*Ny,1)),Nx,Ny,Nf));
%[c1 c2 c3 vec(1:Nx*Ny*Nf)]
choose_flip = repmat(reshape(abs(c2 - c3) < abs(c1 - c3),Nx,Ny,Nf),[1 1 1 Nc]);
tie = repmat(reshape((abs(c2 - c3) == abs(c1 - c3)) & (c1 ~= c2)  ,Nx,Ny,Nf),[1 1 1 Nc]); 
flip_weights = choose_flip + 0.5*tie;
Yfull = sfull_data(:,:,:,:,1).*(1-flip_weights) + sfull_data(:,:,:,:,2).*flip_weights;

% Vectorize data
function vx = vec(x)
vx = x(:);

% Embed data
function y = embed(x,mask)

% Parse inputs
if isempty(mask)
    % No mask
	y = x;
    return;
end
if ~islogical(mask)
    % Uses logical indexing
    error('embed: mask must be logical');
end

% Parse input data
dimx = size(x);
pL   = prod(dimx(2:end));
cl   = class(x);
if issparse(x)
    cl = 'double';
end

% Initialize y
if islogical(x)
	y = false([numel(mask) pL]);       % [np *L]
else
	y = zeros([numel(mask) pL], cl);   % [np *L]
end

% Embed x
if pL > 1
    pLC = num2cell(pL);
	y(mask(:),:) = reshape(x, [], pLC{:});
else
	y(mask(:),1) = x;
end

% Reshape output
if (ndims(mask) == 2) && (size(mask,2) == 1)
    y = reshape(y, [size(mask,1) dimx(2:end)]);   % [(Nd) (L)]
else
    y = reshape(y, [size(mask) dimx(2:end)]);     % [(Nd) (L)]
end

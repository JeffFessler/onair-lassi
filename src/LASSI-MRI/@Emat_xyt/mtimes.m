function y = mtimes(A,x)

if A.adjoint
    % y = A' * x
    arg0 = fft2nc(bsxfun(@times,x,A.mask));
    arg1 = sum(bsxfun(@times,arg0,conj(A.b1)),4);
    arg2 = squeeze(sum(abs(A.b1).^2,4));
    y    = bsxfun(@rdivide,arg1,arg2);
    if ~isempty(A.siz)
        % Reshape into matrix
        y = reshape(y,[],A.siz(end));
    end
else
    % y = A * x
    if ~isempty(A.siz)
        % Reshape into native size
        x = reshape(x,A.siz);
    end
    y = bsxfun(@times,ifft2nc(bsxfun(@times,x,A.b1)),A.mask);
end

% 2-D FFT (normalized, centered)
function FX = fft2nc(X)
[m,n,~] = size(X);
FX = (1 / sqrt(m * n)) * fftshift(fft2(fftshift(X)));

% 2-D IFFT (normalized, centered)
function X = ifft2nc(FX)
[m,n,~] = size(FX);
X = sqrt(m * n) * fftshift(ifft2(fftshift(FX)));

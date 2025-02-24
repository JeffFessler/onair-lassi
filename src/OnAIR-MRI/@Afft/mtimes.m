function y = mtimes(A,x)

if A.adjoint
    % y = A' * x
    y = ifft2n(x .* A.mask);
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
    y = fft2n(x) .* A.mask;
end

% 2-D FFT (normalized)
function FX = fft2n(X)
[m,n,~] = size(X);
FX = (1 / sqrt(m * n)) * fft2(X);

% 2-D IFFT (normalized)
function X = ifft2n(FX)
[m,n,~] = size(FX);
X = sqrt(m * n) * ifft2(FX);

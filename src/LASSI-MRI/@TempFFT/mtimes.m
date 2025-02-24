function y = mtimes(T,x)

if T.adjoint
    % y = T' * x
    y = ifft(ifftshift(x,T.dim),[],T.dim) * sqrt(size(x,T.dim));
    if ~isempty(T.siz)
        % Reshape into matrix
        y = reshape(y,[],T.siz(end));
    end
else
    % Y = T * X
    if ~isempty(T.siz)
        % Reshape into native size
        x = reshape(x,T.siz);
    end
    y = fftshift(fft(x,[],T.dim),T.dim) / sqrt(size(x,T.dim)); 
end

function err = computeErrorMetrics(X,Xtrue)
% Syntax: err = computeErrorMetrics(X,Xtrue);

% Parse inputs
nt = size(Xtrue,3);

% Functions
NRMSEfcn = @(X,Xtrue) norm(X(:) - Xtrue(:)) / norm(Xtrue(:));
pSNRfcn  = @(X,Xtrue) 10 * log10(max(abs(Xtrue(:)))^2 / ...
                                 mean(abs(X(:) - Xtrue(:)).^2));

% Aggregate metrics
NRMSE = NRMSEfcn(X,Xtrue);
pSNR = pSNRfcn(X,Xtrue);

% Per-frame metrics
pNRMSE = nan(1,nt);
ppSNR = nan(1,nt);
for k = 1:nt
    pNRMSE(k) = NRMSEfcn(X(:,:,k),Xtrue(:,:,k));
    ppSNR(k) = pSNRfcn(X(:,:,k),Xtrue(:,:,k));
end

% Return metrics
err.NRMSE  = NRMSE;
err.pSNR   = pSNR;
err.pNRMSE = pNRMSE;
err.ppSNR  = ppSNR;

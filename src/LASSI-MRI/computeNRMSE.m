function nrmse = computeNRMSE(Xhat,Xtrue,ROI)
% Syntax:   nrmse = computeNRMSE(Xhat,Xtrue);
%           nrmse = computeNRMSE(Xhat,Xtrue,ROI);

% Parse inputs
if exist('ROI','var') && ~isempty(ROI)
    % Extract ROI
    Xhat  = Xhat(ROI.rows,ROI.cols,:);
    Xtrue = Xtrue(ROI.rows,ROI.cols,:);
end

% Compute NRMSE
denom = norm(Xtrue(:));
if isnan(denom)
    nrmse = nan;
elseif denom == 0
    nrmse = 0;
else
    nrmse = norm(Xhat(:) - Xtrue(:)) / denom;
end

function Xroi = addROI(X,ROI,lim)
% Syntax:   Xroi = addROI(X,ROI,lim);

% Constants
BORDER = [0, 224, 0];

% Convert to uint8
Xroi = max(0,abs(X) - lim(1));
Xroi = min(1,Xroi / diff(lim));
Xroi = repmat(im2uint8(Xroi), [1, 1, 3]);

% Add ROI border
for i = 1:3
    Xroi(ROI.rows(1)  ,ROI.cols,i) = BORDER(i);
    Xroi(ROI.rows(end),ROI.cols,i) = BORDER(i);
    Xroi(ROI.rows,ROI.cols(1)  ,i) = BORDER(i);
    Xroi(ROI.rows,ROI.cols(end),i) = BORDER(i);
end

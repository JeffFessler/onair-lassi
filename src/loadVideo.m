function X = loadVideo(path,idim,T,dt)
% Syntax: X = loadVideo(path);
%         X = loadVideo(path,idim);
%         X = loadVideo(path,idim,T,dt);

% Parse inputs
RESIZE_VIDEO = (nargin >= 2);
CLIP_FRAMES = (nargin >= 4);

% Load data
data = load(path);
fields = fieldnames(data);
X = im2double(data.(fields{1}));

% Resize video, if necessary
if RESIZE_VIDEO
    if any(~isnan(idim(1:2)))
        X = imresize(X,idim(1:2));
    end
    X = X(:,:,1:nanmin(size(X,3),idim(3)));
end

% Clip frames, if necessary
if CLIP_FRAMES
    nt = size(X,3);
    X = X(:,:,1:(nt - mod(nt - T,dt)));
end

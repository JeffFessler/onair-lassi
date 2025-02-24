function X0t = synthFrameInit(Yslice,Agen,M,Xhat,t,dt,T)
% Syntax:   X0t = synthFrameInit(Yslice,Agen,M,Xhat,t,dt,T);

% Parse inputs
[ny, nx, ~] = size(Xhat);

% Last frame of previous minibatch
tlast = t + T - 1 - dt;
xlast = Xhat(:,:,tlast);
Alast = Agen(true(ny,nx,1),[ny, nx, 1]);
ylast = Alast * xlast; % syntesized data

% Fill missing k-space with synthesized data
Mnew = M(:,:,(tlast + 1):(t + T - 1));
Ynew = Yslice((tlast + 1):(t + T - 1));
Ylast = repmat(ylast,[1, 1, dt, 1]);
missing = repmat(~Mnew,[1, 1, 1, size(Ynew,4)]);
Ynew(missing) = Ylast(missing);

% Synthesize new frames
Anew = Agen(true(ny, nx, dt),[ny, nx, dt]);
Xnew = reshape(Anew' * Ynew,[ny, nx, dt]);

% Synthesized initialization
X0t = cat(3,Xhat(:,:,t:tlast),Xnew);

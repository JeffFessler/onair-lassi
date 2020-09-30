function X0 = zeroFill(Xhat,Yt,Mt,t,dt)
% Syntax: X0 = zeroFill(Xhat,Yt,Mt,t,dt);
%
% Xhat: online estimation (ny x nx x nt)
%   Yt: current time data (ny x nx x T)
%   Mt: current time observation mask (ny x nx x T)
%    t: current time
%   dt: temporal stride
%

% Parse inputs
T = size(Yt,3);
tLast = t + T - 1 - dt; % Last estimated frame

% Fill zeros in new data
Ynew = Yt(:,:,(T + 1 - dt):T);
Mnew = Mt(:,:,(T + 1 - dt):T);
Xrec = repmat(Xhat(:,:,tLast),[1, 1, dt]);
Ynew(~Mnew) = Xrec(~Mnew);

% Construct estimate
X0 = cat(3,Xhat(:,:,t:tLast),Ynew);

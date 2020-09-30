function X0 = interpFill(Xhat,Xit,t,dt)
% Syntax: X0 = interpFill(Xhat,Xit,t,dt);
%
% Xhat: online estimation (ny x nx x nt)
%  Xit: current time interpolated data (ny x nx x T)
%    t: current time
%   dt: temporal stride
%

% Parse inputs
T = size(Xit,3);
tLast = t + T - 1 - dt; % Last estimated frame

% Construct estimate
X0 = cat(3,Xhat(:,:,t:tLast),Xit(:,:,(T + 1 - dt):T));

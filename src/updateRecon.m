function Xhat = updateRecon(Xhat,Xnew,gamma,t,dt)
% Syntax: Xhat = updateRecon(Xhat,Xnew,gamma,t,dt);
%
%  Xhat: reconstruction (ny x nx x nt)
%  Xnew: new data (ny x nx x T)
% gamma: forgetting factor
%     t: current time
%    dt: time stride
%
% Special cases:
%
%     gamma == 0: latest frame
% 0 < gamma  < 1: gamma-weighted average
%     gamma == 1: unweighted average
%

% Parse inputs
T = size(Xnew,3);

% Update estimates
for i = 1:T
    ti = t + i - 1;
    ci = getCount(ti,t,dt,T);
    if gamma == 1
        alpha = 1 / ci;
        beta = (ci - 1) / ci;
    else
        alpha = (1 - gamma) / (1 - gamma^ci);
        beta  = gamma * (1 - gamma^(ci - 1)) / (1 - gamma^ci);
    end
    Xhat(:,:,ti) = alpha * Xnew(:,:,i) + beta * Xhat(:,:,ti);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of times frame ti has been estimated up to (and including) time t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ci = getCount(ti,t,dt,T)
ci = nnz(bsxfun(@plus,(1:dt:t)',0:(T - 1)) == ti);

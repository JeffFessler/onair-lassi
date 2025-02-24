function [X, cost] = cg(Y,X,A,At,D,Dt,Lambda,Lam1,Lam2,Lam3,Lam4,Wx,Wy,Wt,mu1,mu2,beta1,beta2,nIters)
% Conjugate gradient solution to
% \min_X \|A(X) - b\|^2 + \mu\|R(X)\|_1

cost = [];
lam1 = 0.5 * mu1 * beta1;
lam2 = 0.5 * mu2 * beta2;

LamTV  = {Lam1; Lam2; Lam3};
LamTV1 = Lam1 + 1i * 1e-18;
LamTV2 = Lam2 + 1i * 1e-18;
LamTV3 = Lam3 + 1i * 1e-18;
Lam4   = Lam4 + 1i * 1e-18;
for i = 1:nIters
    resY   = A(X) - Y;    
    resL   = X - Lambda; 
    Dx     = D(X);
    resTV1 = Dx{1} - Wx;
    resTV2 = Dx{2} - Wy;
    resTV3 = Dx{3} - Wt;
    resTV  = {resTV1; resTV2; resTV3};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute cost
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eYY   = norm(resY(:))^2;
    eNN   = lam1 * norm(resL(:))^2 + mu1 * abs(Lam4(:)' * resL(:));
    eTV   = lam2 * (norm(resTV1(:))^4 + norm(resTV2(:))^4 + norm(resTV3(:))^4) + ...
             mu2 * abs(conj(LamTV1(:)' * resTV1(:)) + conj(LamTV2(:)' * resTV2(:)) + conj(LamTV3(:)' * resTV3(:))); 
    cost(end + 1) = eYY + eNN + eTV; %#ok
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Line search
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Conjugate gradient direction
    tmpgn = At(resY) + lam1 * (X - Lambda) + lam2 * Dt(resTV);
    gn    = 2 * tmpgn + mu1 * Lam4 + mu2 * Dt(LamTV);
    
    % Search direction  
    if i == 1
        sn    = gn;                                          
        oldgn = gn;
    else
        sn    = gn + (norm(gn(:)) / norm(oldgn(:)))^2 * sn; 
        oldgn = gn;
    end
    
    Asn = A(sn);  
    Dsn = D(sn);
    
    numer =         Asn(:)' * resY(:)  + ...
            lam1 * ( sn(:)' * resL(:)) + ...
       0.5 * mu1 * ( sn(:)' * Lam4(:));
    for k = 1:3
         numer = numer + lam2 * sum(vec(conj(resTV{k}) .* Dsn{k})) + ...
                    0.5 * mu2 * sum(vec(conj(LamTV{k}) .* Dsn{k})); 
    end
    
    denom = norm(Asn(:))^2 + lam1 * norm(sn(:))^2;
    for k = 1:3 
       denom = denom + lam2 * sum(vec(conj(Dsn{k}) .* Dsn{k})); 
    end
    if denom < 1e-18
        break;
    end
    alpha = -real(numer) / real(denom);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % X update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X = X + alpha * sn;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vectorize data
function x = vec(X)
x = X(:);

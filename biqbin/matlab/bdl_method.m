function [ F, G, X, gamma, fopt, t, bestx] = bdl_method( L, b, A, F, G, X, gamma, fopt, t, bestx, itmax, prnt);
% standard  bundle version
% input: bundle (F,G,X), current point: gamma
% call: [ F, G, X, gamma, fopt, t, bestx] = bdl_method( L, b, A, F, G, X, gamma, fopt, t, bestx, itmax, prnt);

if nargin <=11; prnt = 0; end   % no output
% problem data: 
n = size( L,1);         % problem size
m = length( gamma);         % number of constraints
k = length(F);                  % current bundle size

for cnt = 1:itmax               % main iteration

% compute lam, eta and dgamma 
beta = -F -G'*gamma;
[dgamma, lam, eta, t] = lam_eta( beta, G, gamma, t);

% make a step in the direction of -dgamma
gamma_test = gamma + dgamma;   % step into negative subgrd-dir

% compute function at gamma_test
[ftest, xtest, gtest] = fct_eval( gamma_test, L, b, A);

% compute primal X (as convex combination)
xi_lam = X*lam;

% update f and gamma, if progress is made
lmax = max(lam);
del = beta'*lam +dgamma'*dgamma/(2*t) + gamma'*eta +fopt;

% serious step test
if ftest < fopt - 0.05*del
    I = find(lam> lmax*.001);
    G = [G(:,I) gtest];
    X = [X(:,I) xtest(:)];
    F = [F(I); L(:)'*xtest(:)];
    gamma = gamma_test;
    bestx = xi_lam;
    fopt = ftest;
    t = t*1.001;
else             % null step
    I = find( lam(1:k-1)> lmax*.001);
    G = [G(:,I) gtest G(:,k)];    % gamma als letztes glied 
    X = [X(:,I) xtest(:) X(:,k)];
    F = [F(I); L(:)'*xtest(:); F(k)];
    t = t/1.001;
end

k = length( F); 

% display some info
if prnt ~= 0
if mod(cnt,1)==0 | cnt==1 | cnt==itmax;
    lmax = max(gamma);
    ia = length( find( gamma > .001*lmax));  % active constr
    fprintf(' %3.0d  %12.3f %12.3f %8.5f %7.5f %3.0d %5.0d\n', ...
[cnt  fopt ftest norm(dgamma/t) t k ia]);
end
end

end

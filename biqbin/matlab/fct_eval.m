function [f, x, g] = fct_eval( gamma, L, b, A);
% evaluates max-cut function with A(X) <= b and X in Elliptope
% f = <L -A'(gamma),x>  + b'*gamma    function value
% g = b - A(X)                        subgradient
% call:  [f, x, g] = fct_eval( gamma, L, b, A); 

n = size( L,1);                % problem size
m = length(b);
if m > 0;
  Atg = A'*gamma;    % A^t(gamma)
  Atg = reshape(Atg,n,n);  Atg= .5 * (Atg + Atg');
  L0 = L + Atg;                % cost matrix
else
  L0 = L;
end

 [x, dual] = mc_psdpk( L0, 2, 1);         % solve sdp  (L0, 5.5, 1) 

 dualvalue = sum(dual);
if m > 0;
%  f = x(:)'* L0(:) + b'*gamma;  % function value
  f = dualvalue - b'*gamma;
  ax = A*x(:); % A(X)
  g = -b + ax;                   % subgradient  
else
  %f = x(:)'*L0(:); 
  f = dualvalue;
  g = [];
end

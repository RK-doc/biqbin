function [cost, cut] = mc_gw( L, v, trials, II)
% goemans-williamson rounding
% input:  L ... cost matrix (Laplacian), v ... cholesky factor of 
%                                        primal matrix (v'*v=X)
%        (optional) trials ... number of restarts (default=20)  
% output: a cutvector cut of weight cost
% call [ cost, cut] = mc_gw( L, v, trials, II);

if nargin==2; trials = 20; end;
[ n, n1] = size(L);
cost = - 10^10;
k = size(v,1);

for j = 1:trials              % number of rounds is trials
    r = rand( k,1) - .5;      % random vector
	y = sign(v'*r);           % should be +1 or -1
	I = find(y==0);           % should be empty (but just in case)
	if length(I)>0; 
        li = length(I);y(I) = ones(li,1); end;
%  take rounded solution y and improve it through 1-opt
    [ feas, x] = mc_1opt( L, y, II);   
    if feas > cost,
	     cost = feas; cut = x; 
    end;
end;

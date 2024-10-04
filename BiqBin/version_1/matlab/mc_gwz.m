function [ fh, Xh] = mc_gwz( L, X, xh)
% generate cut
% input:  L ... cost matrix (Laplacian), 
%         X ... primal solution
%         xh ... some cut (e.g. xh=rcut(n));
% output: fh, Xh (some solution of value fh)
% uses generalized goemans-williamson rounding
% call [ fh, Xh] = mc_gwz( L, X, xh);

n = length( xh);
its = ceil( n/2);
for i=1:n
    I{i} = find( L(:,i));   % preprocessing
end
fh = xh'*L*xh; done = 0;
while done<2;    % not yet done
    done = done + 1;  
    const = .3+ .6*rand;
    Xs = (1-const)*X + const*xh*xh';
    v = chol( Xs);
    [f0, cut] = mc_gw(L, v, its, I); 
    if f0 > fh %+ 1e-10;
      fh = f0;  
%      fprintf('  better cut, const: %15.4f  %5.3f\n',fh,const);
      xh = cut;
      %const = .3+ .6*rand; 
      %Xs = (1-const)*X + const*xh*xh';  
%      Xs = .5*X + .5*xh*xh';
      done = 0;     % do another round
    end
end    
Xh = cut;


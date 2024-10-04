function [d, lam, eta, t] = lam_eta( beta, G, gamma, t);
% solve lambda eta problem for fixed beta, G, gamma
%  call: [d, lam, eta] = lam_eta( beta, G, gamma, t);
  
k = length( beta);   % number of subgrads 
m = length( gamma);  % number of constraints
eta = zeros(m,1);    % initial value 
done = 1; 
dir_pre = 0; 
lam_cnt = 0;
Q=t*G'*G;

while (done == 1),
        lam_cnt = lam_cnt + 1;
	c=beta-t*G'*eta;
	lam = solve_lambda(c,Q);
	tmp=G*lam;

	eta = max( 0, -gamma/t + tmp);

	d_curr = t*(eta-tmp);
	dir_curr = norm(d_curr);
	if abs(dir_pre-dir_curr)/(1+dir_curr) < 1e-5;
	  done = 0; 
	end;
	if lam_cnt > 50; done = 0; 
%	  disp(['max iter.  ',num2str( abs(dir_pre-dir_curr) )]); 
	  t = t*.95;
	end
	dir_pre = dir_curr;
      end; % end of while
d = d_curr; 

function [ lam, z, y]=solve_lambda(c,Q);
%
% solves
%  (QP)  min <lam,c> + 1/2 <lam, Q lam>
%        s.t. sum(lam_i) = 1, lam >= 0.
% 
%  with input: c, Q  
% output:  lam...optimal solution point
% dual: 
% (D)   max y - 1/2 <lam, Q lam>
%       s.t. z = c + Q lam - e y	       
%       lam >= 0, z >=0 , y free
%  
% call: [ lam, z, y]=solve_lambda(c,Q);

% size of problem  
  k = size(Q,1);
  e = ones(k,1);

% starting triplet (infeasible)
  lam = e/k;
  Qlam= Q*lam;
  tmp = Qlam + c;
  mintmp = min(tmp);
  if mintmp > 1; 
    y= 0; z = tmp;
  else
    y = mintmp - 1; z = tmp - y;
  end

  digits =10.0;%7.0;		 % desired accuracy (6.0)
  mu = (z' * lam) / k;   % calculate barrier parameter mu
  mu=1/2*mu;
  iter = 0;
 
  A_lam=e'*lam;          % linear operator A_lam
  AT_y=y*e;              % adjoint operator AT_y
  
  Res_p=1-A_lam;         % primal residual
  Res_d=c-AT_y-z+Qlam;   % dual residual
 
  D_cost = y - 1/2*lam'*Qlam;	      % initial dual cost
  P_cost = lam'*c + 1/2 *lam'*Qlam;   % intial primal cost
  pd_gap = P_cost-D_cost;	      % pd_gap represents the duality gap
				      % controls stopping condition
  % 			main loop
  % ****************************************************************

  done = 1;
  while done>0;  % while desired accuracy not attained

  iter = iter+1; 
  M = ones(k+1);
  M(k+1,k+1) = 0;
  M(1:k,1:k) = -Q;
%  M = [ -Q e;
%	 e' 0];

 for i=1:k,
   M(i,i)= M(i,i)-z(i)/lam(i);
 end;
% M(1:k,1:k) = M - diag(z./lam); % zum checken
 
  rhs = zeros(k+1,1);     % right hand side predictor
  rhs(1:k) = c-AT_y+Q*lam;%-mu./lam;
  rhs(k+1) = Res_p;

  dw=M\rhs;
  dlam=dw(1:k);
  dy=dw(k+1);

  AT_y=dy*e;  % adjoint operator O^t(dy)
	
  dz=1./lam .* (mu*e-lam.*z - z.*dlam);  % evaluate dz_p

% find steplengths alpha_p and alpha_d s.t.
% lam+alpha_p*dlam >= 0 and z+alpha_d*dz >= 0

  alpha_p=max(-dlam./lam);
  if alpha_p > 0, alpha_p=min(.99/alpha_p,1); else alpha_p=1;end;
  alpha_d=max(-dz./z);
  if alpha_d > 0, alpha_d=min(.99/alpha_d,1); else alpha_d=1;end;

% update

  lam = lam + alpha_p * dlam;
  y = y + alpha_d * dy;
  z = z + alpha_d * dz;
  
  A_lam=e'*lam;  % linear operator A_lam

  AT_y=y*e;  % adjoint operator AT_y
  
  Res_p=1-A_lam;        % primal residual
  Qlam = Q*lam;
  Res_d=c-AT_y-z+Qlam; % dual residual  
    
  mu = (z' * lam) / k;  % evaluate the new mu
  mu=.4*mu;  % reduce mu
if alpha_p+alpha_d > 1.8; mu = mu*.2; end
  
  D_cost = y - 1/2*lam'*Qlam;        % dual cost
  P_cost = lam'*c + 1/2 *lam'*Qlam;  % primal cost
  pd_gap = P_cost-D_cost;	      % duality gap
  if abs(pd_gap) < max( 1, abs(P_cost)) * 10^(-digits); done=0; end
if iter > 30; done = 0; disp(' max iter in solve_lambda...'); end
  
  % display current iteration
% o1=iter;o2=alpha_p;o3=alpha_d;o4=log10(mu);o6=P_cost; o7=D_cost; 
% fprintf('%3d    %5.4f  %5.4f  %5.4f  %5.4f  %5.4f \n',o1, ...
%		   o2,o3,o4,o6,o7);

  end;                  % end of main loop    







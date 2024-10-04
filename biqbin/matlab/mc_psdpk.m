function [X, y, iter, secs] = mc_psdpk( L, digits, silent, X, y);
% input: L ... symmetric  matrix
% output: X   ... primal matrix
%         y   ... dual variables
% solves: max tr(LX): diag(X)=e, X psd, min e'y: Diag(y)-L=Z psd
% call:   [ X, y, iter, secs] = mc_psdpk( L, digits, silent, X, y);
% first input necessary, rest is optional
% silent <> 0  means no output   
% f. rendl, 4/03

start= cputime;
% initialize data
[n, n1] = size( L);             % problem size
b = ones( n,1 );                

if nargin == 1; digits = 2; silent = 0; end;
if nargin == 2; silent = 0; end;
if nargin < 4;
   X = diag( b);
   y = sum( abs(L) )' + 1.;
else
  [dummy, px]=chol(X);
  [dummy, pz]=chol(diag(y)-L);
  if norm(diag(X)-b)> 1e-5; error('diag(X) bad ...'); end
  if max(px, pz)>0; error('input not psd ...');end 
end

Z = diag( y) - L;
phi = b' * y;
psi = L(:)' * X(:);           
delta = phi-psi;

mu = Z(:)' * X( :) /( 2*n);     
alphap = 1; alphad = 1; iter = 0; cholcnt = 0;

if silent == 0;
disp(['      iter    alphap    alphad    log(gap)  lower      upper ']);
fprintf(' %2.0d   %10.3f %10.3f   %8.5f  %10.3f  %10.3f\n', ...
	 [ iter alphap alphad log10(delta) psi phi]); 

end;

while delta > 10^(-digits) %max([abs(phi) 1]) * 10^(-digits)  % while duality gap too large

        Zi = inv( Z); iter = iter + 1;
        dzi = diag(Zi);
        Zi = (Zi + Zi')/2;
% predictor step solves: Z * X + diag( dy1) * X + Z * dX1 = 0
        dy1 =  (Zi .* X) \ ( - b);       % solve for dy
        tmp = scale_mat( Zi, dy1);
	dX1 = - tmp * X -X;
	dX1 = (dX1 + dX1')/2;
% corrector step
% solves: diag(dy2)*X + Z*dX2-muI + diag(dy1)*dX1 = 0
	rhs = mu *dzi - (Zi .* dX1)*dy1;
	dy2 = (Zi .* X) \ rhs;
	tmp = scale_matr( dy2, X);    % diag( dy2) * X
	tmp1 = scale_matr( dy1, dX1);  % diag(dy1) * dX1
	dX2 = mu*Zi - Zi*( tmp + tmp1);
% final steps
        dy = dy1 + dy2;
        dX = dX1 + dX2;
	dX = ( dX + dX')/2;     % symmetrise

%	dX(1:n+1:end) = zeros(n,1);
	
% find steplengths alphap and alphad
        alphap =  1;
        [Zi,posdef] = chol( X + alphap * dX);
        cholcnt = cholcnt + 1;
        while posdef ~= 0,
                alphap = alphap * .8;
                [Zi,posdef] = chol( X + alphap * dX);
                cholcnt = cholcnt + 1;
                end;
% stay away from boundary
        if alphap < 1, alphap = alphap * .95; end;
%	alphap = alphap*.95;
        X = X + alphap * dX;

        alphad = 1;
        dZ = sparse(diag( dy));
        [Zi,posdef] = chol( Z + alphad * dZ);
        cholcnt = cholcnt + 1;
        while posdef ~= 0;
                alphad = alphad * .8;
                [Zi,posdef] = chol( Z + alphad * dZ);
                cholcnt = cholcnt + 1;
                end;
        if alphad < 1, alphad = alphad * .95; end;
%        alphad = alphad*.95;

% update
        y = y + alphad * dy;
        Z = Z + alphad * dZ;
        mu = X( :)'* Z( :) /(2*n);
        if alphap + alphad > 1.6,
        mu = mu * .5; end;
        if alphap + alphad > 1.9,
        mu = mu/(5); end;   % reduce mu, if stepsize good
        phi = b' * y;
        psi = L( :)'* X( :);         
        delta = phi-psi;

%	if min([alphap alphad]) < 1e-3; disp(' konvergenz??');
%	   end
% display current iteration
if silent == 0;
  fprintf(' %2.0d   %10.3f %10.3f   %8.5f  %10.3f  %10.3f\n', ...
	 [ iter alphap alphad log10(delta) psi phi]); 
end;
        end;            % end of main loop
secs = cputime - start;





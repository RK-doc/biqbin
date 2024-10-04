function [X, y, iter, secs] = mc_psd( L, digits, silent);
% input: L ... symmetric  matrix
% output: X   ... primal matrix
%         y   ... dual variables
% solves: max tr(LX): diag(X)=e, X psd, min e'y: Diag(y)-L=Z psd
% call:   [ X, y, iter, secs] = mc_psd( L, digits, silent);
% first input necessary, rest is optional
% silent <> 0  means no output   
% f. rendl, 2/99

start = cputime;
% initialize data
[n, n1] = size( L);             % problem size
if nargin == 1; digits = 2; silent = 0; end;
if nargin == 2; silent = 0; end;
b = ones( n,1 );                
X = diag( b);
y = sum( abs(L) )' + 1.1;
Z = diag( y) - L;
phi = b' * y;
psi = L(:)' * X(:);           
delta = phi-psi;

mu = Z(:)' * X( :) /( 4*n);     
alphap = 1; alphad = 1; iter = 0;

if silent == 0;
disp([' iter      alphap    alphad   log(gap)    lower      upper ']);
fprintf(' %2.0d   %10.3f %10.3f   %8.5f  %10.3f  %10.3f\n', ...
	 [ iter alphap alphad log10(delta) psi phi]); 
end;

while delta > 10^(-digits) %max([abs(phi) 1]) * 10^(-digits)  % while duality gap too large

        Zi = inv( Z); iter = iter + 1;
        dzi = diag(Zi);
        Zi = (Zi + Zi')/2;
        dy =  (Zi .* X) \ (mu * dzi - b);       % solve for dy
        tmp = zeros(n);        
        for j=1:n; 
	  tmp(:,j) = Zi(:,j)*dy(j); 
	end;
        dX = -tmp * X + mu*Zi -X;
%       dX  = - Zi * diag( dy) * X + mu * Zi - X;
        dX = ( dX + dX')/2;     % symmetrise

% find steplengths alphap and alphad
        alphap =  1;
        [Zi,posdef] = chol( X + alphap * dX);
        while posdef ~= 0,
                alphap = alphap * .8;
                [Zi,posdef] = chol( X + alphap * dX);
                end;
% stay away from boundary
        if alphap < 1, alphap = alphap * .95; end;
        X = X + alphap * dX;

        alphad = 1;
%        dZ = sparse(diag( dy));
        dZ = diag( dy);      % for octave
        [Zi,posdef] = chol( Z + alphad * dZ);
        while posdef ~= 0;
                alphad = alphad * .8;
                [Zi,posdef] = chol( Z + alphad * dZ);
                end;
        if alphad < 1, alphad = alphad * .95; end;

% update
        y = y + alphad * dy;
        Z = Z + alphad * dZ;
        mu = X( :)'* Z( :) /(2*n);
        if min(alphap,alphad) < .5; mu = mu * 1.5; end
        if alphap + alphad > 1.6,
        mu = mu * .75; end;
        if alphap + alphad > 1.9,
        mu = mu/(1.+.1 * iter); end;   % reduce mu, if stepsize good
        phi = b' * y;
        psi = L( :)'* X( :);         
        delta = phi-psi;

% display current iteration
if silent == 0;
      fprintf(' %2.0d   %10.3f %10.3f   %8.5f  %10.3f  %10.3f\n', ...
	 [ iter alphap alphad log10(delta) psi phi]); 
end;
        end;            % end of main loop
secs = cputime - start;






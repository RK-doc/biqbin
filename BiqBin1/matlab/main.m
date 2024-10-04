function [bnd,bestcut] = main( L, var);
% input: L      .. Laplacian
% call: [ gamma, bX,bestcut, A, bnd, X, secs]= main( L, var);

% compute first bound and initial set of triangs and 
% intialize bundle

%set variant to 3 (=triangles only) if only one input
if nargin==1; var=3; end
tstart = tic; n = size( L,1);
[b,T, F,G,X,A,gamma,fopt,t,bestx,bestcut]=mcstart( L);
cutval = bestcut'*L*bestcut;
bnd = fopt;
if abs( bnd-cutval<.99)
    return;
end

maxv = 1; maxv5 = 1; maxv7 = 1; % to initialize stopping condition

% define some default setting
% for ising istances: maxit = 200, bdl_it = 5; cutting planes = 100;
maxit= 20;%15;     % number of iterations   
bdl_it = 10;   % iterations of the bundle method for fixed set of constraints
level = 0.01;   % violation level to stop early 

% number of bundle iterations
it = bdl_it; f5(1) = 0; f5old = 0;
%fprintf('        time       bnd     cut  bdl_it  bdl    m     viol3   nc3   viol5   nc5 \n');

% main loop
for cnt = 1:maxit

% increase bundle iterations after the first rounds
%if cnt==11; it = it + 2; end

% call standard bundle
[F,G,X,gamma,fopt,t,bestx]=bdl_method(L,b,A,F,G,X,gamma,fopt,t,bestx, it);
bnd = fopt;  % keep best bound
bX = reshape( bestx, n, n);   % bestx as matrix

% purge inactive constraints
[ b, A, T, G, gamma] = bdl_purge( b, A, T, G, gamma);

% call heuristic with new bestx
startcut = rcut(n);
[fh, xh] = mc_gwz( L, bX, startcut);
if fh > cutval;
   bestcut = xh; cutval = fh; 
end

now = toc(tstart); 
bsize = length( F);   % bundle size
m = length(b);
%fprintf('%3.0d %8.2f %10.3f %10.3f %3.0d %4.0d %5.0d * ',cnt,now, fopt, cutval,it,bsize,m); 

% stopping condition
if var == 3; done = maxv<level;  end
if var == 5; done = maxv5<level; end
if abs( fopt-cutval<.95); done= 1>0; end   % we assume integer data, so stop if gap smaller than 1

if cnt == maxit | done,  fprintf( '\n'); secs=toc(tstart); return, end

% separate new triangles
[ A3new, T, gammanew] = separation( bX, T);

% new data
A = [A; A3new];
b = ones( size(A,1),1);
tnew = length(gammanew); 
maxv = max(gammanew);   % largest triangle violation
if length(gammanew)<1; maxv=0.09; tnew=0; end
%fprintf( '%6.3f %5.0d', maxv,tnew);

if maxv<.4 & var>3;       % .15 % .2
   [C5,f5] = sep_c5( bX);
   [A5] = C5_to_a( n, C5');
   gamma5 = 1- A5*bX(:);
   maxv5 = max( gamma5);   % largest C5 violation
   %fprintf(' %7.3f %5.0d ', 1-f5(1), length(f5));
%   fprintf('\n');
   A = [A; A5];
   b = ones( size(A,1), 1);
   gammanew = [gammanew; gamma5];
else
%   fprintf('\n');
end
%fprintf('\n');   

   

% make new gamma
gamma = [gamma; 0.005*t*gammanew];   % 0.1 might be too big

% evaluate at new point
[ f, x, g] = fct_eval( gamma, L, b, A);

% update bundle info (F,G)
ng = length(g);     % number of constraints
r = size( X, 2);    % bundle size 
G = zeros( ng, r+1);   % initialize subgradients
for i=1:r
   G(:,i) = -b + A*X(:,i);
end

% now store old center and add new point
G = [ G(:,1:r-1) g, G(:,r)];
F = [ F(1:r-1); L(:)'*x(:); F(r)];
X = [ X(:,1:r-1) x(:) X(:,r)];

% new estimate for t
%t = 0.3 * (fopt - cutval) / (g'*g);
t = 1.05*t;     % 1.1*t
secs = toc(tstart);
end

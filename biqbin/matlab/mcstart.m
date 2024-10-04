function [ b, T, F, G, X, A, gamma, fopt, t, bestx, bestcut] = mcstart( L);
% setup for max-cut and bundle with triangle inequalities
% (F,G,X) bundle, t ... bundle parameter
% call:  [ b, T, F, G, X, gamma, fopt, t, bestx, bestcut] = mcstart( L);

% first solve basic sdp
Xs = mc_psd(L, 3, 1); 
f = Xs(:)'*L(:);

% generate first cut
n = size( L,1);
xh = rcut( n);    % a random cut
[fh, xh] = mc_gwz( L, Xs, xh);


% first results
%fprintf(' start:   bound: %12.3f     cut: %12.3f \n', f, fh);

% check pruning condition
fopt = f;
if abs( fopt-fh<.99)
    b = [];
    T = [];
    F = [];
    G = [];
    X = Xs;
    bestx = X;
    A = [];
    gamma = [];
    t = 0;
    bestcut = xh;
    return;
end

% separate triangles
[ A, T, g] = separation( Xs,[]);
b = ones( size(A,1), 1);
%fprintf(' triangles added: %5.0d,  violation: %7.3f %7.3f \n', length(b), min(g), max(g));

% starting t
t = 0.5 * (f - fh) / (g'*g);  % 0.3 

% starting gamma
gamma = 0.1 * t * g;

% first evaluation at gamma
[ f, x, g] = fct_eval( gamma, L, b, A);
%fprintf( ' starting bound: %12.3f \n', f);

% setup for bundle
F = L(:)'*x(:); G = g; X = x(:); 
fopt = f; 
bestx=X; bestcut = xh;

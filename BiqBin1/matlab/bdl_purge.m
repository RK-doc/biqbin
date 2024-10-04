function [ b, A, T, G, gamma]  = bdl_purge( b, A, T, G, gamma)
% purge inactive constraints
% call:  [ b, A, T, G, gamma]  = bdl_purge( b, A, T, G, gamma);

% determine inactive constraints, given gamma
maxg = max(gamma);
mtri = size( T,1);
I = find( gamma > 0.0001*maxg);      % dual variables should be large enough
Itri = find(gamma(1: size(T,1)) > .0001*maxg);
% fprintf(' purge: %5.0d triangs kept from %5.0d \n', length( I), length( gamma));

gamma = gamma( I);
b = b( I);
A = A( I, :);
T = T( Itri, :);
G = G( I, :);

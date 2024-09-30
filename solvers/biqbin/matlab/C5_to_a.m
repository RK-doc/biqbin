function [ A] = C5_to_a( n, C5);
% call [ A] = C5_to_a( n, C5);

m = size( C5,1);  % number of sets
A = sparse( [],[],[],n*n,m, m*25);
e = ones(1,5);

for i=1:m
   raw = C5( i,1:5);
   sup1 = abs( raw);
   el1 = sign( raw);
   I = kron(e, sup1);
   J = kron( sup1, e);
   S = kron(el1,el1);
   A0 = sparse( I, J, S, n,n, 25);
   A(:,i) = A0(:);
end
A = A';

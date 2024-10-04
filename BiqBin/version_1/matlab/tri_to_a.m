function [A,b] = tri_to_a( n, sup, el);
% call [ A,b ] = tri_to_a( n, sup, el);

m = size( sup,1);  % number of triangle inequ
A = sparse( [],[],[],n*n,m, m*9);
e = ones(1,3);

for i=1:m
   sup1 = sup( i,:);
   el1 = el( i,:);
   I = kron(e, sup1);
   J = kron( sup1, e);
   S = kron(el1,el1);
   A0 = sparse( I, J, S, n,n, 9);
   A(:,i) = A0(:);
end
A = A';
b=ones(m,1);

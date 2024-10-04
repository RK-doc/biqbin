function Y = scale_matr(y, X);
% compute Y = diag(y) * X
% y is vector, X is matrix, i.e. row scaling of X by y
% call:   X = scale_matr( y, X);

n = length( y);  
n1 = size( X,2);
Y = zeros(n, n1);   
for j=1:n; 
  Y(j,:) = X(j,:)*y(j); 
end;
 
function X = scale_mat(Z, y);
% compute X = Z * diag(y)
% y is vector, Z is matrix, i.e. column scaling of Z by y
% call:   X = scale_mat( Z, y);

n = length( y);  
X = zeros(n);        
for j=1:n; 
  X(:,j) = Z(:,j)*y(j); 
end;
  
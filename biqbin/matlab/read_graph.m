function [L,A] = read_graph(file)
% CONVERT reads aa graph in edge list format from file 
% and produces Laplacian matrix L 

e = dlmread(file);
n = e(1,1);
m = e(1,2);
e = e(2:end,:);

% data for max-cut problem
B = sparse(e(:,1),e(:,2),e(:,3),n,n);
A = B + B' - diag(B);

L = 1/4*(diag(A*ones(n,1))-A);

end


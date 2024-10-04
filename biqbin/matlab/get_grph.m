function [e,w] = get_grph(A,filename,ID)
% Given the (sym) matrix A it outputs
% the weighted  graph with edges e and weights w.
% if A(i,j) = 0, then no edge is generated
% 
% call:  [ e, w] = get_grphInt(A);
  
[n,n1] = size(A);
B=abs(A)>0;
m = sum(sum(B + diag(diag(B))))/2;
e = zeros(m,2);
w = zeros(m,1);
k = 0;
for i = 1:n;
for j = i:n;
         if A(i,j) ~= 0;
                k = k+1;
                e(k,:) = [i j];
                w(k) = A(i,j);
                end;
end; end;

A = [e , w];

instance = sprintf('./data/%s.%s',ID,filename);
fileID = fopen(instance,'wt');
fprintf(fileID,'%d %d\n',n,m);
fprintf(fileID,'%d %d %.f \n',A'); % do I need more decimals %.6f?
fclose(fileID);

end

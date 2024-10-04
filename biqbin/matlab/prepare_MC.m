function prepare_MC(rpath,filename)
% MATLAB function for transformation phase of BiqBin, i.e.
%
% min x'Fx + c'x  s.t.  Ax = b,  x in {0,1}^n
% 
% to max-cut instance

% a) read general BQP with linear constraints
instance = sprintf('%s%s',rpath,filename);
[A,b,c,F] = biqbin2matlab(instance);

% b) compute penalty parameter and obtain max-cut instance
[val,upp,feas] = matlab2graphMIN(A,b,c,F,filename);

% c) save val (offset), upp (upper bound to check infeasibility of original
% problem as val - maxcut > upp) and feas \in {-1,0,1} (1 feasible, -1
% original problem infeasible, 0 compute optimal cut and check condition 
% val - maxcut > upp for infeasibility

% save offset
fid = fopen('./data/data.txt','w');
fprintf(fid,'feasible\n%d\n',feas);
fprintf(fid,'offset\n%f\n',val);
fprintf(fid,'upp\n%f\n',upp);
fclose(fid); 


end


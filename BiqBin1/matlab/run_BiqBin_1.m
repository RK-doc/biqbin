% MATLAB script for transformation phase of BiqBin, i.e.
%
% min x'Fx + c'x  s.t.  Ax = b,  x in {0,1}^n
% 
% to max-cut instance

% a) read general BQP with linear constraints
[A,b,c,F] = biqbin2matlab('example_feasible.txt'); %input file as function input!!


% %% test
% % DATA
% n = 20;             % for big n cplex takes long
% 
% % number of constraints
% m = 10;
% 
% c = randi([-1,1],n,1);
% 
% F = randi([-1,1],n,n);
% F = F'*F;           
% 
% A = randi([-1,1],m,n);
% 
% % let x be the feasible solution solution
% x = randi(2,n,1);
% x = mod(x,2);           % 0/1 vector
% 
% %rhs
% b = A*x;

%%
% b) compute penalty parameter and obtain max-cut instance
[val,upp,feas] = matlab2graphMIN(A,b,c,F);

% read graph from file
L = read_graph('MC_instance.txt');

% c) run BiqBin (in this case main)
[~, cut]= main(L,5);
maxcut = cut'*L*cut;

% d) transform the solution back to 0-1 variables
% by deleting the last component in using linear transformation 
cut = cut(1:end-1);
cut = (cut + ones(size(cut)))/2;

%the optimum solution (note if biqbin did not finish you can still output
%approximate solution)
if feas == 1 % original problem is surely feasible
    fprintf('optimal value: %f\noptimal cut:\n', val - maxcut);
    %disp(cut');
else
    if val - maxcut > upp
        fprintf('the problem is infeasible\n');
    else
        fprintf('optimal value: %f\noptimal cut:\n', val - maxcut);
        %disp(cut');
    end
end


% obtain optimal solution of the original problem
% integer programs are solved with cplexmiqp.

% all variables are binary = 'B'
myString = '';
n = length(c);

for i=1:n
    myString = strcat(myString,'B');
end

[xConstrained, valueConstrained,exitflag] = cplexmiqp(2*F,c,[],[],A,b,[],[],[],[],[],myString);
valueConstrained
exitflag



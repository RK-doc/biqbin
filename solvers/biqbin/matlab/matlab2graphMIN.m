function [val,upp,feas] = matlab2graphMIN(A,b,c,F,filename)

% Given the parameters of the problem with all entries integer, it 
% creates graph with integer weights. It outputs the value and 
% (possibly) the bounds, the parameters and the Laplacian matrix.
%
% A,b,c,F       input parameters (0/1 formulation)
% val           offset value
% L             Laplacian matrix
% pen           penalty paraemter
% time          time
% low and upp   bounds used to compute pen

%% Setting of the problem data and {-1,1} transformation
tic; 
[m,n] = size(A); 
b1 = b - (A*ones(n,1))/2;   A1 = A/2;
c1 = (c+(ones(1,n)*F)')/2;  F1 = F/4;
alpha = sum(c/2) + sum(sum(F/4));

% Definition of matrices (of size n+1) for objective and penalized const.
OBJ = [ alpha , c1'/2 ; c1/2 , F1 ];                % objective part
A2 = A1'*A1;    b2 = b1'*b1;    lin = -A1'*b1;
CON = [ b2 , lin' ; lin , A2 ];                     % constraint part

%% Upper Bound (Goemans-Williamson heuristic, otherwise MET \cap X_I)
penTriv = 2*sum(sum(abs(F)));       % a (trivial) valid penalty parameter
Ql = OBJ + penTriv*CON;             % +because min
Dl = diag(sum(Ql)); 
GW = Dl - Ql; 
Xs = mc_psd(GW,5,1); cut = rcut(n+1); 
[bnd,cut] = mc_gwz(GW,Xs,cut);      % Goemans-Williamson heuristic

if cut'*CON*cut==0                  % condition for feasibility
    ub = -bnd + sum(sum(Dl));       % upper bound if feasible instance
    feas = 1;
else %0 is a feasible vector
    M = [b1, -A1];                    N = null(M);
    [nc,nv] = size(N);              % n+1 constraints, n+1-rk(M) variables
    FF = N'*[ alpha , (c1/2)' ; (c1/2) , F1 ]*N; % Objective
    for i = 1:nc
        NN{i} = N(i,:)'*N(i,:); % constraints
    end
    env = ones(1,nv); % all one vector (e) of size nv
    enc = ones(1,nc); % e of size nc
    clear prob
    prob.bardim = nv;
    % Objective
    Fl = tril(FF); % lower triangular part of matrix FF
    [ro,co,vo] = find(Fl);
    enzF = ones(1,length(ro)); % e of size #nonzero el. of tril(FF)
    prob.barc.subj = enzF;
    prob.barc.subk = ro';
    prob.barc.subl = co';
    prob.barc.val  = vo';
    % Constraints   
    prob.bara.subi = [];
    prob.bara.subj = [];
    prob.bara.subk = [];
    prob.bara.subl = [];
    prob.bara.val  = [];
    prob.a         = sparse([],[],[],nc,nc);
    for i=1:nc
        aus = tril(NN{i});
        [raus,caus,vaus] = find(aus);
        enzaus = ones(1,length(raus)); % e of size #nonzero el. of tril(NNi)
        prob.bara.subi = [prob.bara.subi,i*enzaus];
        prob.bara.subj = [prob.bara.subj,enzaus];
        prob.bara.subk = [prob.bara.subk,raus'];
        prob.bara.subl = [prob.bara.subl,caus'];
        prob.bara.val  = [prob.bara.val ,vaus'];
    end
    % Bounds of constraints
    prob.blc = enc;
    prob.buc = enc;
    
    % no output of mosek
    param.MSK_IPAR_INFEAS_REPORT_LEVEL = 0;
    param.MSK_IPAR_LOG_INFEAS_ANA = 0;
    param.MSK_IPAR_LOG = 0;
    param.MSK_IPAR_INFEAS_REPORT_AUTO = 'MSK_OFF';
    
    
    % Optimization part and results
    [~,res] = mosekopt('maximize',prob,param);
    sol = res.sol.itr;
    ub  = sol.pobjval;
    
    % check feasiblity of SDP
    if strcmp(sol.prosta,'PRIMAL_INFEASIBLE')
       feas = -1;  
       val = 0;
       upp = 0;
       return;
    else
       feas = 0;
    end
end

%% Lower bound (with 5-clique ineq, i.e., code of Franz)
    AUS = diag(sum(OBJ));   
	MAT = OBJ - AUS;
	bnd = main(-MAT,5); % 3 for triangles
	lb = -bnd + sum(sum(AUS));
    
%% Fixing the bounds (integer) and defining \sigma
pen = floor(ub-lb) + 1;                     % integer, pen = u-l+\eps ok
OBJt = [ F1, c1/2 ; c1'/2 , alpha ];        % objective part, for Timotej
CONt = [ A2 , lin ; lin' , b2 ];            % constraint part, for Timotej
low = floor(lb); upp = ceil(ub);            % integer bounds
L = OBJt + pen*CONt;                        % Q in Lasserre paper                   
val = sum(sum(L));                          % offset
BM = L - diag(diag(L)); 
get_grph(4*BM,filename);  % create graph, 4 makes integer
%fid = fopen('val.txt','w');	fprintf(fid,'%.2f',val);	fclose(fid); % save offset
time = toc;
end
%fid = fopen('feas.txt','w'); fprintf(fid,'%d',feas); fclose(fid); % save
%feasibility
%
% SOLUTION IS: val - cut
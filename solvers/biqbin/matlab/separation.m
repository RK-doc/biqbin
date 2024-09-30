function [ A, T, gamma_new] = separation(x, Told)
% x is matrix
% generate new inequalities 
% call:  [ A, T, gamma_new] = separation( x, Told);
    
  n = size(x,1);

% number of new inequalities is new_ineq
  new_ineq = n*10; % war 4n fuer sdp  % n*n/10
 
% find new violated constraints at x
  X1 = (ones(n)-x)/2;         % transform into 0-1 model
  y = reshape( X1, n^2, 1);   % y vector for X
 
  [Tn, g_new] = trianglesep( y, n, new_ineq);
% [Tn, dummy1, hashn, g_new] = tri_sep2( y, n, new_ineq, h_tmp);

  m = length( g_new);         % number of new ineq.
  if m == 0, 
    disp('nothing new violated.'); 
    gamma_new = []; return    % return with old data structure
  end;
  if m > 0;
    Tn = reshape( Tn, m, 4);
   end

% merge new with old constraints
  mold = size(Told,1);  % old number of constraints
  T = [Told; Tn]; mtot= size(T,1); 
  [T1, IA] = unique( T, 'rows'); 
   Ikeep = intersect( mold+1:mtot, IA);
  Ikeep = Ikeep-mold;      % these are now the new triangles
  m = length( Ikeep);
  Tn = Tn(Ikeep,:); 
  T = [Told; Tn];
  
 % generate A(X) <=b
 [sup,el] = to_ip( Tn);
 A = tri_to_a( n, sup, el);
 gamma_new = 1 - A*x(:);

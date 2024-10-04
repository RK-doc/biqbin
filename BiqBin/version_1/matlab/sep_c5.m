function [I5, f] = sep_c5( X);
% call: [I5, f] = sep_c5( X);
% use qap

e = ones(5,1);H1 = e*e';
e(1) = -1; H2 = e*e';
e(2) = -1; H3 = e*e';
trials = 300;
f = zeros( 3*trials, 1); 
I5 = zeros(5, 3*trials); 
for i=1:trials
   [perm,val] = qap_simul2_c( H1, X);
   f(i) = val; 
   I5(:,i) = perm(1:5); 
   [perm,val] = qap_simul2_c( H2, X);
   f(trials + i) = val; 
   perm(1) = -perm(1);
   I5(:,trials + i) = perm(1:5); 
   [perm,val] = qap_simul2_c( H3, X);
   f(2*trials + i) = val;  
   perm(1) = -perm(1); perm(2) = -perm(2); 
   I5(:,2*trials + i) = perm(1:5); 
end
% first check for identical sets
I5keep = I5;
I5 = sort( abs(I5));   % sort columnwise
[dummy, K] = unique(I5', 'rows');
I5 = I5keep(:, K);
f = f(K);
% now remove nonviolations and sort 
K = find(f<0.99); f = f(K); I5=I5(:,K);
[dummy, K] = sort(f);
f = f(K);
I5 = I5(:,K);




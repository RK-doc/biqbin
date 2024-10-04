function [cost, xnew] = mc_1opt( L, x, II);
% improves a given cut x with Laplace matrix L (cost=x'*L*x)
% new cut is xnew
% call: [cost, xnew] = mc_1opt( L, x. II)

Lx = L*x;       % auxiliary vector
d = diag(L);
cost = x'*Lx;
delta = d - x.*Lx;
[best, i] = max( delta);
while best > 0.00001
        I = find(L(:,i));
         %I = II{i};
%        Lx(I) = Lx(I)  - 2 * x(i)*L(I,i);     % update L*x;
        if x(i)>0;
           Lx(I) = Lx(I)  - 2 *L(I,i);
                else
           Lx(I) = Lx(I)  + 2 *L(I,i);
                end;
%        Lx = Lx  - 2 * x(i)*L(:,i);     % update L*x;
        x( i) = -x( i);                 % update new cut
        cost = cost + 4*best;           % update weight of cut
        delta = d - x.*Lx;              % update new differences
        [best, i] = max( delta);        % find new champion
        end;
xnew = x;

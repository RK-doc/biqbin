function x = rcut(n)
% generate a random cut x ( i.e. +-1 vector of size n
% call x = rcut( n)

x = ones(n,1);
for i=1:n
        y= rand;
        if y >= .5; x( i) = -1; end
end;
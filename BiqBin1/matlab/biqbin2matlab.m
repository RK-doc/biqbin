function [A,b,c,F] = biqbin2matlab(fname)

fid = fopen(fname);

% n = number of variables:
line = fgetl(fid); 
aus = str2num(line);
n = aus(1);
m = aus(2);

% initialization of matrices
A = zeros(m,n);
b = zeros(m,1);
c = zeros(n,1);
F = zeros(n,n);


line = fgetl(fid);
% defition of A
line = fgetl(fid);
while line(1)~='b'
    aa = str2num(line);
    A(aa(1),aa(2)) = aa(3);
    line = fgetl(fid);
end

% definition of b
line = fgetl(fid);
while line(1)~='F' & line(1)~=' F'
    bb = str2num(line);
    b(bb(1),1) = bb(2);
    line = fgetl(fid);
end

% definition of F
line = fgetl(fid);
while line(1) ~= 'c'
    FF = str2num(line);
    F(FF(1),FF(2)) = FF(3);
    line = fgetl(fid);
end
F = F + F' - diag(diag(F));

% definition of c
line = fgetl(fid);
while length(line) > 1
    cc = str2num(line);
    c(cc(1),1) = cc(2);
    line = fgetl(fid);
end

fclose(fid);

end

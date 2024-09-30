function [sup, el] = to_ip( T_g);
% wechsel zu interior point daten struktur
% call :  [sup, el] = to_ip( T_g);
  
m = size( T_g,1);   % gruber's datenstruktur

sup = T_g( 1:m, 1:3);
el = ones( m, 3);
for i = 1:m;
  if T_g(i,4) == 1; el(i, 3) = -1; end
  if T_g(i,4) == 2; el(i, 2) = -1; end
  if T_g(i,4) == 3; el(i, 1) = -1; end
end  
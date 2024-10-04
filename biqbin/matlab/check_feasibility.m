function [ output_args ] = check_feasibility(filename)
% checks before and anfter maxcut computations if the original binary QP
% with linear constraints is feasible or not

fid = fopen(filename);

% second line of data contains number indicating feasbility
fgetl(fid);
feas = fscanf(fid,'%d');

% 3 options: 
% feas = 1 --> original problem is feasible


end


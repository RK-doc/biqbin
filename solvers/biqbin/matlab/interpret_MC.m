function interpret_MC(rpath,filename)
% rpath = path to the folder containing filename
% filename = output of max-cut solver in JSON format

% read data.txt
dataFile = fopen('./data/data.txt');
fgetl(dataFile);
feas = fscanf(dataFile,'%d');
fgetl(dataFile);
offset = fscanf(dataFile,'%f');
fgetl(dataFile);
upp = fscanf(dataFile,'%d');


% change name of outputfile
filename = sprintf('%s%s',rpath,filename);

% read results from BiqBin maxcut solver in JSON file
% and check (in)feasibility of original problem
fid = fopen(filename);

temp_sol = sprintf('%s%s',rpath,'temp_sol');
foutput = fopen(temp_sol,'w');

% read lines and print them to foutput
line = fgetl(fid);fprintf(foutput,line);fprintf(foutput,'\n');
line = fgetl(fid);fprintf(foutput,line);fprintf(foutput,'\n');
line = fgetl(fid);fprintf(foutput,line);fprintf(foutput,'\n');
line = fgetl(fid);fprintf(foutput,line);fprintf(foutput,'\n');

% read 5th line for optimum value
line = fgetl(fid);
line = strtrim(line);
maxcut = sscanf(line,'"Solution": %d');

fprintf(foutput, '\t\t\"Value\":');

if feas == 1 % original problem is surely feasible
    fprintf(foutput, ' %d,\n', offset - maxcut);
else
    if offset - maxcut > upp
		fclose(foutput);
		foutput = fopen(temp_sol,'w');
        fprintf(foutput,'{"Message": "The problem is infeasible."}');
    else
        fprintf(foutput, ' %d,\n', offset - maxcut);
    end
end

% read 6th line for optimum cut
line = fgetl(fid);
line = strtrim(line);
ex = extractBetween(line,': ( ', ' )');
cut_index = str2num(cell2mat(ex));

% read 13th line for number of original vertices + 1
fgetl(fid); % 7
fgetl(fid); % 8
fgetl(fid); % 9
fgetl(fid); % 10
fgetl(fid); % 11
fgetl(fid); % 12
line = fgetl(fid);
line = strtrim(line);
n = sscanf(line,'"Vertices": %d');

% start reading fid from beginning
fclose(fid);

fid = fopen(filename);
fgetl(fid); % 1
fgetl(fid); % 2
fgetl(fid); % 3
fgetl(fid); % 4
fgetl(fid); % 5
fgetl(fid); % 6

% transform the solution back to 0-1 of original problem:
% last vertex is deleted and indeces from max-cut solver are set to 0, so 1 -> 0 and 0 -> 1
cut = ones(n-1,1);
for i = 1:length(cut_index)
	cut(cut_index(i)) = 0;
end

fprintf(foutput, '\t\t\"Solution\": \"[ ');
for i = 1:n-1
	fprintf(foutput, '%d ', cut(i));
end
fprintf(foutput, ']\",\n');

% copy some more lines
line = fgetl(fid);fprintf(foutput,line);fprintf(foutput,'\n');
line = fgetl(fid); % root node bound is skipped and not written to output file
line = fgetl(fid);fprintf(foutput,line);fprintf(foutput,'\n');
line = fgetl(fid);fprintf(foutput,line);fprintf(foutput,'\n');
line = fgetl(fid);fprintf(foutput,line);fprintf(foutput,'\n');

% ... and finish
fprintf(foutput, '\t\t\"ObjectTypeString\": \"BQP\",\n');
fprintf(foutput,'\t\t\"NumVariables\": %d\n', n-1);
fprintf(foutput, '\t}\n');
fprintf(foutput, '}\n');

fclose(fid);
fclose(foutput);


end


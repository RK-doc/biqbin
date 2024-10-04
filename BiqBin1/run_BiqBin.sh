#!/bin/bash

unset LD_PRELOAD
module load MATLAB/2019a

biqbin=/home/hrgat/Namizje/General_BiqBin_parallel_argv_final/
user_folder=/home/hrgat/Namizje/General_BiqBin_parallel_argv_final/
rpath=$biqbin/Instances/
filename=example.txt
ID=TEST

# remove output file if already exists
if [ -f $user_folder$ID.output ] ; then
	rm $user_folder$ID.output
fi

# run matlab script to prepare max-cut instance
matlab -nodisplay -nosplash -nodesktop -r "addpath('/home/hrgat/mosek/8/toolbox/r2014a');cd('$biqbin./matlab/');prepare_MC('$rpath','$filename');exit"


# check for infeasibility
line=$(sed -n '2p' $biqbin./matlab/data/data.txt)

if [[ $line == -1 ]]
then
	echo '{"Message": "The problem is infeasible."}' > $ID.output 
else

	# run biqbin max-cut solver on the same file name located in matlab/data
	sbatch run_maxcut.sh $filename $ID

	# check when max-cut solver finished!
	until [ -f $biqbin$ID.output ]
	do
		sleep 5
	done

	# remove temporaty files from data folder (data.txt must NOT be deleted)
	rm $biqbin./matlab/data/$filename
	
	# prepare output file with transformed results
	matlab -nodisplay -nosplash -nodesktop -r "cd('$biqbin./matlab/');interpret_MC('$user_folder','$ID.output');exit"

	# delete max-cut output and rename temp file
	rm $userfoled$ID.output
	
	#name=temp_sol
	mv $user_folder./temp_sol $user_folder$ID.output
	
fi






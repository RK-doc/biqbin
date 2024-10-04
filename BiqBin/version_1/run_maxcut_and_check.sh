#!/bin/bash

biqbin=$1
user_folder=$2
filename=$3
ID=$4

unset LD_PRELOAD
source ${biqbin}setupenv.sh
module load MATLAB/2019a

# check for infeasibility
line=$(sed -n '2p' ${biqbin}matlab/data/$ID.data.txt)

echo 'Line: '$line

if [[ $line == -1 ]]; then
	echo '{"Message": "The problem is infeasible."}' > ${user_folder}$ID.output 
else

	# run biqbin max-cut solver on the same file name located in matlab/data

	# sbatch $biqbin/run_mc_web.sh $ID.$filename $ID
	mpirun -n 12 ${biqbin}biqbin ${biqbin}matlab/data/$ID.$filename ${biqbin}params ${user_folder}$ID 5 30 1	

	# check when max-cut solver finished!
	until [ -f ${user_folder}$ID.output.tmp ]
	do
		sleep 5
	done

	# remove temporary files from data folder
	rm ${biqbin}matlab/data/$ID.$filename
	
	# prepare output file with transformed results
	matlab -nodisplay -nosplash -nodesktop -r "cd('${biqbin}matlab/');interpret_MC('$user_folder','$ID');exit"

	# delete max-cut output and rename temp file
	rm ${user_folder}$ID.output.tmp
	
fi

# delete data file from matlab folder
rm ${biqbin}matlab/data/$ID.data.txt

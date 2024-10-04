#!/bin/bash

biqbin=$1
user_folder=$2
filename=$3
ID=$4

# Load libraries
unset LD_PRELOAD
module load MATLAB/2019a

# remove output file if already exists
if [ -f ${user_folder}$ID.output ] ; then
	rm ${user_folder}$ID.output
fi

# run matlab script to prepare max-cut instance
matlab -nodisplay -nosplash -nodesktop -r "addpath('~/mosek/8/toolbox/r2014a');cd('${biqbin}matlab/');prepare_MC('$user_folder','$filename','$ID');exit"





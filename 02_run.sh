#!/bin/bash

# Use this script to submit jobs on cluster that uses schedular like SLURM 

# If you are running locally, just run job.sh scripts in each run directory. 

totruns=5
run=1
while [ $run -le $totruns ];do
    cd ${run}
    if [ $1 ]; then
	if [ $1 == "sbatch" ]; then
	    sbatch job.sh
	elif [ $1 == "qsub" ]; then
	    qsub job.sh
	fi
    else
	printf "Running: run - $run \n"
	bash job.sh > job.out
	wait
    fi
    cd ..
    run=$((run+1))
done
echo "Done!"
exit

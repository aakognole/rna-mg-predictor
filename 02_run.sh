#!/bin/bash

# Use this script to submit jobs on cluster that uses schedular like SLURM 

# If you are running locally, just run job.sh scripts in each run directory. 

totruns=5
run=1
while [ $run -le $totruns ];do
    cd ${run}
    if [ $1 ]; then
	sbatch job.sh
    else
	qsub job.sh
    fi
    cd ..
    run=$((run+1))
done

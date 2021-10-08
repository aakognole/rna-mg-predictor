#!/bin/bash
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

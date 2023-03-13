#!/bin/bash

NumForSim=11 #number of jobs for each simulation - 1

qsub ./FDPS-basalt-sphere-collision.sh > job_id.txt
sleep 0.1s
	
for ((i=0 ; i<$NumForSim ; i++))
do
    JobIDSTR=$(<./job_id.txt)
    JobID=${JobIDSTR%.*}
    echo $JobID
    qsub -W depend=afterok:$JobID ./mid-script.sh > job_id.txt
    sleep 0.1s
done

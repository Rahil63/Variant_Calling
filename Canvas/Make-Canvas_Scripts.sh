#!/bin/bash

## Script produces executable scripts for each BAM pair in $(pwd) to be run with Canvas,
## submitted to SLURM:

for i in *-t_*.bam; do 
  BASENAME=${i%-t_SortedRmdup.bam}
  cat slurm_header.txt \
    <(echo '#SBATCH --job-name=Canvas_'${BASENAME}) \
    <(echo '#SBATCH --output='${BASENAME}'_canvas.log') \
    <(echo 'BASENAME='${BASENAME}) \
    Canvas_master.sh \
    > Canvas_run_${BASENAME}.sh; done

##
chmod +x Canvas_run_*.sh  

##
for i in Canvas_run_*.sh
  do
  sbatch $i
  done

#!/bin/bash

## Loop through Manta or Strelka directories and submit the jobs to SLURM,
## where start_manta.sh is simply the following lines:

#######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=70
#SBATCH --partition=normal
#SBATCH --time=08:00:00 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a_toen03@uni-muenster.de
#SBATCH --job-name=manta_sv
#SBATCH --output=manta.log
#######


for i in MantaDir_RG*
  do 
  cp start_manta.sh 
  ./$i/Somatic
  cd ./$i/Somatic
  sbatch start_manta.sh
  cd ../../
  cp start_manta.sh
  ./$i/Germline
  done
  
exit

for i in StrelkaDir_RG*
  do 
  cp start_strelka.sh 
  ./$i/Somatic
  cd ./$i/Somatic
  sbatch start_strelka.sh
  cd ../../
  cp start_strelka.sh
  ./$i/Germline
  done  
  

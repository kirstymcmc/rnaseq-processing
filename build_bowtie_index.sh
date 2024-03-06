#!/bin/bash

#SBATCH --qos bbdefault
#SBATCH --account plackarg-spl-bioinformatic-analysis
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --ntasks 30 # request 130 cores for the job.
#SBATCH --time 2-00:00:00 # this requests 2 days
#SBATCH --mail-type ALL 


#Load required modules 
module purge 
module load bear-apps/2022b
module load Bowtie2/2.5.1-GCC-12.2.0

cd ../data/test_dir/SILVA

bowtie2-build SILVA.db silva --threads 30 --large-index
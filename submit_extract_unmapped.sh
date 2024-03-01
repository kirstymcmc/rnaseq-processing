#!/bin/bash

#SBATCH --qos bbdefault
#SBATCH --account plackarg-spl-bioinformatic-analysis
#SBATCH --ntasks 12 # request 12 cores for the job. N.B. check whether the fastq-dump can parallelise, else this is redundant and you should set to "1"
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 1000 # this requests 2 hours, but you will need to adjust depending on runtime. Test job execution time with just a couple of input files then scale accordingly
#SBATCH --mail-type ALL 


module purge;
module load bluebear


#run from inside full_dataset/data folder


../rnaseq-processing/extract_unmapped.sh . ./2_trimmed ./4_unmapped "u"


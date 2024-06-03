#!/bin/bash
#SBATCH --qos bbdefault
#SBATCH --account plackarg-spl-bioinformatic-analysis
#SBATCH --ntasks 30 # request 30 cores for the job. N.B. check whether the fastq-dump can parallelise, else this is redundant and you should set to "1"
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 2-0 # this requests 2 days
#SBATCH --mail-type ALL 

set -e

module purge;
module load bluebear



# Move to directory containing clean reads 
cd ../data/6_rrna_filtered
# Ensure trinity image is contained within 6_rrna_filtered folder

apptainer exec -e trinityrnaseq.v2.15.1.simg  Trinity \
          --seqType fq \
          --samples_file ../trinity_sample_table.txt \
          --max_memory 24G --CPU $SLURM_NTASKS \
          --output ../trinity_out_dir
          
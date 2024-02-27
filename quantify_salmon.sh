#!/bin/bash

#SBATCH --qos castles
#SBATCH --account plackarg-spl-bioinformatic-analysis
#SBATCH --ntasks 8 # request 8 cores for the job. N.B. check whether the fastq-dump can parallelise, else this is redundant and you should set to "1"
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 1000 # this requests 2 hours, but you will need to adjust depending on runtime. Test job execution time with just a couple of input files then scale accordingly
#SBATCH --mail-type ALL 

set -e

module purge;
module load bluebear
module load bear-apps/2021b
module load Salmon/1.9.0-GCC-11.2.0


REF_INDEX=ref/salmon_index
FILES=../data/2_trimmed/*1_trimmo.fq.gz
MAPPED_FILES=data/quants


for f in $FILES
do

f=${f##*/}
f=${f%_1_trimmo.fq.gz}


salmon quant -i $REF_INDEX -l A \
-1 ../data/2_trimmed/${f}_1_trimmo.fq.gz \
-2 ../data/2_trimmed/${f}_2_trimmo.fq.gz \
-p 12 --validateMappings -o $MAPPED_FILES/${f}_quant

done

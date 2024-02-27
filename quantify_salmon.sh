#!/bin/bash

#SBATCH --qos bbdefault
#SBATCH --account plackarg-spl-bioinformatic-analysis
#SBATCH --mem-per-cpu=6750M
#SBATCH --ntasks 54 # request 54 cores for the job. 
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 12:00:00 # this requests 12 hours, 
#SBATCH --mail-type ALL

set -e

module purge;
module load bluebear
module load bear-apps/2022b
module load Salmon/1.10.1-GCC-12.2.0

mkdir -p ../data/3_quants

REF_INDEX=../ref/salmon_index
FILES=../data/2_trimmed/*1_trimmo.fq.gz
MAPPED_FILES=data/3_quants


for f in $FILES
do

f=${f##*/}
f=${f%_1_trimmo.fq.gz}


salmon quant -i $REF_INDEX -l A \
-1 ../data/2_trimmed/${f}_1_trimmo.fq.gz \
-2 ../data/2_trimmed/${f}_2_trimmo.fq.gz \
-p 54 --validateMappings -o $MAPPED_FILES/${f} --writeUnmappedNames \
--gcBias --seqBias

done

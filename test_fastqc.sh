#!/bin/bash

#SBATCH --qos bbdefault
#SBATCH --account plackarg-spl-bioinformatic-analysis
#SBATCH --ntasks 12 # request 12 cores for the job.
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 1-00:00 # this requests 1 day
#SBATCH --mail-type ALL 

module purge;
module load bluebear
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-foss-2019b-Python-3.7.4


# Go to work folder
workdir="../data/test_dir"
#cd $workdir

# Initial quality check using FastQC and MultiQC
echo "Preprocessing pipeline started" > $workdir/pipeline_log.txt
echo Current Date and Time is: `date +"%Y-%m-%d %T"` >> $workdir/pipeline_log.txt

SECONDS=0
echo "Running FastQC on raw reads" | tee -a $workdir/pipeline_log.txt

# run FastQC on all files in the dir trimmed_reads
find $workdir/4_unmapped/ -name '*.fastq' | xargs fastqc

# combine reports with MultiQC in the folder '4_unmapped'
multiqc $workdir/4_unmapped/ -o $workdir/4_unmapped/ -n 4_unmapped_report

echo "Time needed to finish FastQC step on raw reads: $SECONDS seconds" >> $workdir/pipeline_log.txt

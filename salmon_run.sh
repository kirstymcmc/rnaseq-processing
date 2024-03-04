#!/bin/bash

#SBATCH --qos bbdefault
#SBATCH --account plackarg-spl-bioinformatic-analysis
#SBATCH --ntasks 30 # request 30 cores - necessary for parallel step
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 5-00:00 # this requests 5 days
#SBATCH --mail-type ALL 

module purge;
module load bluebear
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-foss-2019b-Python-3.7.4

##################################################
# Initial quality check using FastQC and MultiQC

# Define work folder
workdir="../data/"
# Define the output directory for FastQC reports
fastqc_output_dir="$workdir/qc_raw"
multiqc_output_dir="$workdir/qc_raw"

# Check if the output directory exists, if not, create it
if [ ! -d "$fastqc_output_dir" ]; then
  mkdir -p "$fastqc_output_dir"
fi


echo "Preprocessing pipeline started" > $workdir/pipeline_log.txt
echo Current Date and Time is: `date +"%Y-%m-%d %T"` >> $workdir/pipeline_log.txt

SECONDS=0
echo "Running FastQC on raw reads" | tee -a $workdir/pipeline_log.txt

# run FastQC on all files in the dir trimmed_reads
find $workdir/1_raw/ -name '*.fq.gz' | xargs fastqc --outdir $fastqc_output_dir

# combine reports with MultiQC in the folder '4_unmapped'
multiqc $fastqc_output_dir -o $multiqc_output_dir -n raw_qc_report

echo "Time needed to finish FastQC step on raw reads: $SECONDS seconds" >> $workdir/pipeline_log.txt

##################################################
# Trim raw reads with TrimGalore
# Adapter removal and read quality trimming of paired-read fastq-files

# load required module
module purge;
module load bluebear
module load bear-apps/2019b/live
module load Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4

# create dir for trimmed reads
workdir="../data"
mkdir -p "$workdir/2_trimmed"

# run on all pair-read fastq.gz files in a folder
# number of core fixed to 8! Read cutadapt manual for explanation

fqdir=$workdir/1_raw
for r1 in "$fqdir"/*1.fq.gz; do
    r2=${r1%1.fq.gz}2.fq.gz
    if [[ -f $r2 ]] ; then
trim_galore --paired --retain_unpaired --phred33 \
--output_dir $workdir/2_trimmed --length 50 -q 5 \
--stringency 1 -e 0.1 --cores 8 "$r1" "$r2"
    else
        echo "$r2 not found" >&2
    fi
done

echo "Time needed to finish Trim Galore step: $SECONDS seconds" >> $workdir/pipeline_log.txt

#######################################################
# Re-ren fastqc and multiqc on trimmed reads to assess trimming quality

module purge;
module load bluebear
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-foss-2019b-Python-3.7.4


# Define work folder
workdir="../data/"
# Define the output directory for FastQC reports
fastqc_output_dir="$workdir/qc_trimmed"
multiqc_output_dir="$workdir/qc_trimmed"

# Check if the output directory exists, if not, create it
if [ ! -d "$fastqc_output_dir" ]; then
  mkdir -p "$fastqc_output_dir"
fi



echo Current Date and Time is: `date +"%Y-%m-%d %T"` >> $workdir/pipeline_log.txt
SECONDS=0
echo "Running FastQC on trimmed reads" | tee -a $workdir/pipeline_log.txt

# run FastQC on all files in the dir trimmed_reads
find $workdir/2_trimmed/ -name '*.fq.gz' | xargs fastqc --outdir $fastqc_output_dir

# combine reports with MultiQC in the folder '4_unmapped'
multiqc $fastqc_output_dir -o $multiqc_output_dir -n trimmed_qc_report

echo "Time needed to finish FastQC step on trimmed reads: $SECONDS seconds" >> $workdir/pipeline_log.txt

##################################################
# READS are now ready to quantify with Salmon
#Salmon requires large amount of mem - submit as separate job
# Ensure you have made a salmon index with index_salmon.sh first
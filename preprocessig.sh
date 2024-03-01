#!/bin/bash

#SBATCH --qos bbdefault
#SBATCH --account plackarg-spl-bioinformatic-analysis
#SBATCH --ntasks 12 # request 12 cores for the job.
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 1-00:00 # this requests 1 day
#SBATCH --mail-type ALL 

module purge
module load bluebear
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-foss-2019b-Python-3.7.4


##########################################################################################
# script adapted from 
# https://github.com/matevzl533/Noccaea_praecox_transcriptome/blob/main/preprocessing.sh #
# raw reads must be in folder '4_unmapped' in the $workdir folder                        #
# processed reads will be in the 5_filtered folder                                       # 
##########################################################################################

# Go to work folder
workdir = "../data/test_dir"
cd $workdir

# Initial quality check using FastQC and MultiQC
echo "Preprocessing pipeline started" > $workdir/pipeline_log.txt
echo Current Date and Time is: `date +"%Y-%m-%d %T"` >> $workdir/pipeline_log.txt

SECONDS=0
echo "Running FastQC on raw reads" | tee -a $workdir/pipeline_log.txt

# run FastQC on all files in the dir trimmed_reads
find $workdir/4_unmapped/ -name '*.fastq.gz' | xargs fastqc

# combine reports with MultiQC in the folder '4_unmapped'
multiqc $workdir/4_unmapped/ -o $workdir/4_unmapped/ -n 4_unmapped_report

echo "Time needed to finish FastQC step on raw reads: $SECONDS seconds" >> $workdir/pipeline_log.txt

module purge
module load bear-apps/2021b
module load Python/3.9.6-GCCcore-11.2.0
module load Rcorrector/1.0.5-GCC-11.2.0

# Removing erroneous k-mers from Illumina paired-end reads
# Run with the highest possible number of cores
SECONDS=0
echo "Running Rcorrector on raw reads" | tee -a $workdir/pipeline_log.txt
mkdir cor_reads
echo "Corrected reads will be in cor_reads folder" | tee -a $workdir/pipeline_log.txt

# Raw reads are in folder '4_unmapped'
fqdir=$workdir/4_unmapped
for r1 in "$fqdir"/*1_unmapped.fastq.gz; do
    r2=${r1%1_unmapped.fastq.gz}2_unmapped.fastq.gz
    if [[ -f $r2 ]] ; then
        run_rcorrector.pl -1 $r1 -2 $r2 -t 24 -od '$workdir/cor_reads/'
    else
        echo "$r2 not found" >&2
    fi
done
echo "Time needed to finish Rcorrector step on raw reads: $SECONDS seconds" >> $workdir/pipeline_log.txt

# Discard read pairs for which one of the reads is deemed unfixable
# Original python script from https://github.com/harvardinformatics/TranscriptomeAssemblyTools
# Translated to Python3
# Not memory intensive, therefore assign as many cores as possible with -P => parallel processes

SECONDS=0
echo "Running Filter Uncorrectable on corrected reads" | tee -a $workdir/pipeline_log.txt
echo "Filtered reads will be in folder clean_reads" | tee -a $workdir/pipeline_log.txt
mkdir clean_reads
cd clean_reads
ls $workdir/cor_reads/*1_unmapped.cor.fq.gz | xargs -P30 -I@ bash -c 'python $workdir/FilterUncorrectabledPEfastq.py -1 "$1" -2 ${1%1.*.*}2_unmapped.cor.fq.gz' _ @
echo "Time needed to finish Rcorrector step on raw reads: $SECONDS seconds" >> $workdir/pipeline_log.txt
cd

SECONDS=0
echo "Running Trim Galore" | tee -a pipeline_log.txt
echo Current Date and Time is: `date +"%Y-%m-%d %T"` >> $workdir/pipeline_log.txt

# adapter removal and read quality trimming of paired-read fastq-files
# trim_galore --paired --retain_unpaired --phred33 --output_dir trimmed_reads --length 50 -q 5 --stringency 1 -e 0.1 --cores 2 SAMPLE_R1.fastq.gz SAMPLE_R2.fastq.gz

# load required module
module purge 
module load Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4

# create dir for trimmed reads
mkdir trimmed_reads

# run on all pair-read fastq.gz files in a folder
# number of core fixed to 8! Read cutadapt manual for explanation

fqdir=$workdir/clean_reads
for r1 in "$fqdir"/*1_unmapped.cor.fq; do
    r2=${r1%1.cor.fq}2_unmapped.cor.fq
    if [[ -f $r2 ]] ; then
trim_galore --paired --retain_unpaired --phred33 \
--output_dir $workdir/trimmed_reads --length 50 -q 5 \
--stringency 1 -e 0.1 --cores 8 "$r1" "$r2"
    else
        echo "$r2 not found" >&2
    fi
done

echo "Time needed to finish Trim Galor step: $SECONDS seconds" >> $workdir/pipeline_log.txt

# Check against SILVA rRNA db
# The files you want are *_blacklist_paired_unaligned.fq.gz

#Load required modules 
module purge 
module load bear-apps/2022b
module load Bowtie2/2.5.1-GCC-12.2.0

fqdir=$workdir/trimmed_reads
for r1 in "$fqdir"/*1_unmapped.cor_val_1.fq; do
    r2=${r1%1.cor_val_1.fq}2_unmapped.cor_val_2.fq
    if [[ -f $r2 ]] ; then
        bowtie2 --quiet --very-sensitive-local \
--phred33  -x $workdir/../SILVA/SILVA.db -1 "$r1" -2 "$r2" --threads 6 \
--met-file ${r1%.fq}_bowtie2_metrics.txt \
--al-conc-gz ${r1%.fq}_blacklist_paired_aligned.fq.gz \
--un-conc-gz ${r1%.fq}_blacklist_paired_unaligned.fq.gz  \
--al-gz ${r1%.fq}_blacklist_unpaired_aligned.fq.gz \
--un-gz ${r1%.fq}_blacklist_unpaired_unaligned.fq.gz
    else
        echo "$r2 not found" >&2
    fi
done

SECONDS=0
echo "Running FastQC on processed reads" | tee -a $workdir/pipeline_log.txt
echo Current Date and Time is: `date +"%Y-%m-%d %T"` >> $workdir/pipeline_log.txt

# Moving the files we want to 5_filtered folder

mkdir $workdir/5_filtered
mv $workdir/trimmed_reads/*_paired_unaligned.fq* $workdir/5_filtered
cd $workdir/5_filtered

# Postprocessing quality check
# run FastQC on all files in the trimmed_reads folder
#load required modules 

module purge
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-foss-2019b-Python-3.7.4

find -name '*_paired_unaligned.fq*.gz' | xargs fastqc

# Move Trim Galore and FastQC reports to endQC folder for MultiQC
mv $workdir/trimmed_reads/*.txt $workdir/endQC
mv $workdir/5_filtered/*.html $workdir/endQC

# combine reports with MultiQC
multiqc $workdir/endQC/ -o $workdir/endQC/ -n 5_filtered_report

echo "Time needed to finish second FastQC step on processed reads: $SECONDS seconds" | tee -a $workdir/pipeline_log.txt
echo "FastQC on processed reads finished." | tee -a $workdir/pipeline_log.txt
echo "###################################" | tee -a $workdir/pipeline_log.txt
echo "Preprocessing pipeline finished." | tee -a $workdir/pipeline_log.txt
cd
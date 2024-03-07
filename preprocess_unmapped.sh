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


# Define work folder
workdir="../data"
# Define the output directory for FastQC reports
fastqc_output_dir="$workdir/4_qc"
multiqc_output_dir="$workdir/4_qc"

# Check if the output directory exists, if not, create it
if [ ! -d "$fastqc_output_dir" ]; then
  mkdir -p "$fastqc_output_dir"
fi

##################################################
# Initial quality check using FastQC and MultiQC
echo "Preprocessing pipeline started" > $workdir/pipeline_log.txt
echo Current Date and Time is: `date +"%Y-%m-%d %T"` >> $workdir/pipeline_log.txt

SECONDS=0
echo "Running FastQC on raw reads" | tee -a $workdir/pipeline_log.txt

# run FastQC on all files in the dir trimmed_reads
find $workdir/4_unmapped_63/ -name '*.fq.gz' | xargs fastqc --outdir $fastqc_output_dir

# combine reports with MultiQC in the folder '4_unmapped_63'
multiqc $fastqc_output_dir -o $multiqc_output_dir -n 4_unmapped_report

echo "Time needed to finish FastQC step on raw reads: $SECONDS seconds" >> $workdir/pipeline_log.txt

#######################################################
module purge;
module load bluebear
module load bear-apps/2021b
module load Python/3.9.6-GCCcore-11.2.0
module load Rcorrector/1.0.5-GCC-11.2.0

# Removing erroneous k-mers from Illumina paired-end reads
# Run with the highest possible number of cores
SECONDS=0
echo "Running Rcorrector on raw reads" | tee -a "$workdir/pipeline_log.txt"
mkdir -p "$workdir/5_corrected"
echo "Corrected reads will be in 5_corrected folder" | tee -a "$workdir/pipeline_log.txt"

# Raw reads are in folder '4_unmapped_63'
fqdir="$workdir/4_unmapped_63"
for r1 in "$fqdir"/*1_unmapped.fq.gz; do
    r2="${r1%1_unmapped.fq.gz}2_unmapped.fq.gz"
    if [[ -f $r2 ]] ; then
        run_rcorrector.pl -1 "$r1" -2 "$r2" -t 24 -od "$workdir/5_corrected"
    else
        echo "$r2 not found" >&2
    fi
done
echo "Time needed to finish Rcorrector step on raw reads: $SECONDS seconds" >> "$workdir/pipeline_log.txt"

# Discard read pairs for which one of the reads is deemed unfixable
# Original python script from https://github.com/harvardinformatics/TranscriptomeAssemblyTools
#script updated to python 3 - renamed FilterUncorrectabledPEfastq_p3.py
# Not memory intensive, therefore assign as many cores as possible with -P => parallel processes

module purge;
module load bluebear
module load bear-apps/2021b
module load Python/2.7.18-GCCcore-11.2.0-bare

SECONDS=0
echo "Running Filter Uncorrectable on corrected reads" | tee -a $workdir/pipeline_log.txt
echo "Filtered reads will be in folder 5_cor_filtered" | tee -a $workdir/pipeline_log.txt
mkdir -p "$workdir/5_cor_filtered"

ls $workdir/5_corrected/*1_unmapped.cor.fq | xargs -P30 -I@ bash -c 'python FilterUncorrectabledPEfastq_P3.py -1 "$1" -2 "${1%1_unmapped.cor.fq}2_unmapped.cor.fq"' _ @

#fixed reads will output to folder that Python script is in
# move to desired folder (5_cor_filtered) and remove "unfixrm" prefix

echo "Moving filtered reads to 5_cor_filtered and removing unfixrm_ prefix" | tee -a $workdir/pipeline_log.txt

# Loop through all files starting with 'unfixrm_' in the current directory
for file in unfixrm_*; do
    # Remove the 'unfixrm_' prefix and prepare the new file path
    newfile="${file#unfixrm_}" # This syntax removes the 'unfixrm_' prefix from the filename
    mv "$file" "$workdir/5_cor_filtered/$newfile" #moves new file to 5_cor_filtered
done

echo "Time needed to finish Rcorrector step on raw reads: $SECONDS seconds" >> $workdir/pipeline_log.txt

#######################################################
# Remove rrna contamination with Bowtie2

module purge;
module load bluebear
module load bear-apps/2022b
module load Bowtie2/2.5.1-GCC-12.2.0

echo "Removing rRNA contamination with Bowtie2" | tee -a $workdir/pipeline_log.txt
echo "Filtered reads will be in folder 6_rrna_filtered" | tee -a $workdir/pipeline_log.txt
SECONDS=0

fqdir=$workdir/5_cor_filtered
outdir=$workdir/6_rrna_filtered
mkdir -p $workdir/6_rrna_filtered

for r1 in "$fqdir"/*1_unmapped.cor.fq; do
    r2=${r1%1_unmapped.cor.fq}2_unmapped.cor.fq
        base=$(basename "$r1" | cut -d '_' -f 1-3)

        bowtie2 --quiet --very-sensitive-local \
        --phred33  -x "$workdir/SILVA/silva" -1 "$r1" -2 "$r2" --threads 6 \
        -S "$fqdir/${base}_alignment.sam" \
        --met-file ${r1%.fq}_bowtie2_metrics.txt \
        --un-conc-gz "$outdir/${base}_clean_%.fq.gz"  \

done

echo "Time needed to finish Bowtie2 step: $SECONDS seconds" >> $workdir/pipeline_log.txt
echo Current Date and Time is: `date +"%Y-%m-%d %T"` >> $workdir/pipeline_log.txt

#######################################################
# Final quality check using FastQC and MultiQC
SECONDS=0
echo "Running FastQC on processed reads" | tee -a $workdir/pipeline_log.txt
echo Current Date and Time is: `date +"%Y-%m-%d %T"` >> $workdir/pipeline_log.txt


module purge
module load bluebear
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-foss-2019b-Python-3.7.4

# Define the output directory for FastQC reports
fastqc_output_dir="$workdir/endQC"
multiqc_output_dir="$workdir/endQC"

# Check if the output directory exists, if not, create it
if [ ! -d "$fastqc_output_dir" ]; then
  mkdir -p "$fastqc_output_dir"
fi

find $workdir/6_rrna_filtered -name '*.fq.gz' | xargs fastqc --outdir $fastqc_output_dir

# combine reports with MultiQC
multiqc $workdir/endQC/ -o $workdir/endQC/ -n 7_filtered_report

echo "Time needed to finish second FastQC step on processed reads: $SECONDS seconds" | tee -a $workdir/pipeline_log.txt
echo "FastQC on processed reads finished." | tee -a $workdir/pipeline_log.txt
echo "###################################" | tee -a $workdir/pipeline_log.txt
echo "Preprocessing pipeline finished." | tee -a $workdir/pipeline_log.txt
cd
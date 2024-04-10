#!/bin/bash
#SBATCH --qos bbdefault
#SBATCH --account plackarg-spl-bioinformatic-analysis
#SBATCH --ntasks 8 # request 8 cores for the job. N.B. check whether the fastq-dump can parallelise, else this is redundant and you should set to "1"
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 1000 # this requests 2 hours, but you will need to adjust depending on runtime. Test job execution time with just a couple of input files then scale accordingly
#SBATCH --mail-type ALL 

set -e

module purge;
module load bluebear
module load bear-apps/2022b
module load SAMtools/1.17-GCC-12.2.0
module load HISAT2/2.2.1-gompi-2022b
module load BEDTools/2.30.0-GCC-12.2.0
module load Biopython/1.81-foss-2022b




# Map decoy sequences to the genome using HISAT2 and convert to bed format
#REF_INDEX=../../../hisat2_stuff/reference/hisat_index/index
#FILES=../../data/4_decoys/*_1_unmapped.fq.gz
#MAPPED_FILES=../../data/4_decoys_mapped
#FILEPATH=../../data/4_decoys
#mkdir -p $MAPPED_FILES

#for f in $FILES
#do
#    f=${f##*/}
#    f=${f%_1_unmapped.fq.gz}

#    hisat2 -p 8 \
#    -x $REF_INDEX \
#    --summary-file $MAPPED_FILES/${f}_summary.txt \
#    -2 $FILEPATH/${f}_2_unmapped.fq.gz | \
#    -1 $FILEPATH/${f}_1_unmapped.fq.gz \
#    samtools sort -o $MAPPED_FILES/${f}_sorted.bam

    # Convert sorted BAM to BED
#    bedtools bamtobed -i $MAPPED_FILES/${f}_sorted.bam > $MAPPED_FILES/${f}.bed
#done

# Make exon, intron and intergenic bed files from Crichardii gff
python gff_tools.py

#
bedtools coverage -a exons.bed -b ../../data/4_decoys_mapped/Apx_130d_BR1.bed -f 0.9 -F 0.9 -hist
bedtools coverage -a intergenic.bed -b ../../data/4_decoys_mapped/Apx_130d_BR1.bed -f 0.9 -F 0.9 -hist
bedtools coverage -a introns.bed -b ../../data/4_decoys_mapped/Apx_130d_BR1.bed -f 0.9 -F 0.9 -hist
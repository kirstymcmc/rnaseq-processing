#!/bin/bash

#SBATCH --qos bbdefault
#SBATCH --account plackarg-spl-bioinformatic-analysis
#SBATCH --ntasks 30 # request 30 cores - necessary for parallel step
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 1-00:00 # this requests 1 day
#SBATCH --mail-type ALL 

# Map clean reads to filtered trinity assembly with Bowtie2

module purge;
module load bluebear
module load bear-apps/2022b
module load DIAMOND/2.1.8-GCC-12.2.0

# BLAST trinity contigs against nt database using DIAMOND
# To retrieve taxonomy information (needed to classify alignments by taxon), taxon info
# must be supplied in the form of a taxon map file. This file is a tab-delimited file
# These can be downloaded from NCBI's taxonomy database
# Files are in ref/blastdb/

workdir=".."
##make custom nr database with taxon ids 
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip

#diamond makedb --in nr.gz -d nr --taxonmap prot.accession2taxid.gz --taxonnodes taxdmp.zip --threads 10

diamond blastx --db $workdir/ref/blastdb/nr.dmnd --query $workdir/data/Trinity_cdhit90.fasta \
 --out $workdir/Trinity_cdhit90_nr.blastx \
 --outfmt 6 qseqid sseqid qstart qend sstart send evalue pident mismatch staxids sscinames sskingdoms skingdoms sphylums \
 --max-target-seqs 1 --evalue 1e-10 --threads 30

echo "Time needed to finish Diamond analysis: $SECONDS seconds" >> $workdir/pipeline_log.txt
echo Current Date and Time is: `date +"%Y-%m-%d %T"` >> $workdir/pipeline_log.txt
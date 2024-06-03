#!/bin/bash

#SBATCH --qos bbdefault
#SBATCH --account plackarg-spl-bioinformatic-analysis
#SBATCH --ntasks 30 # request 30 cores - necessary for parallel step
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 1-00:00 # this requests 1 day
#SBATCH --mail-type ALL 
set -e
# Map clean reads to filtered trinity assembly with Bowtie2

# purge and load modules

#module purge;
#module load bluebear
#module load bear-apps/2022b
#module load DIAMOND/2.1.8-GCC-12.2.0
#module load TransDecoder/5.7.1-GCC-12.2.0

#  ________________________________________________________________________________________
# | BLAST trinity contigs against nt database using DIAMOND                                |
# | To retrieve taxonomy information (needed to classify alignments by taxon), taxon info  |
# | must be supplied in the form of a taxon map file. This file is a tab-delimited file    |
# | These can be downloaded from NCBI's taxonomy database                                  |
# | Files are in ref/blastdb/                                                              |
#  ________________________________________________________________________________________

workdir="/rds/projects/p/plackarg-spl-bioinformatic-analysis/full_dataset"


##make custom nr database with taxon ids - only need to run once
#cd $workdir/ref/blastdb
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
#unzip taxdmp.zip
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz

#diamond makedb --in nr.gz -d nr --taxonmap prot.accession2taxid.gz --taxonnodes nodes.dmp --taxonnames names.dmp --threads 30


cd $workdir
mkdir -p ./data/7_blast

#filter trinity contigs to remove seqeunces less than 500bp in length
#perl ./rnaseq-processing/filter_fasta.pl ./data/Trinity_cdhit90.fasta 500 > ./data/Trinity_cdhit90_500.fasta

# Predict ORFs 
#cd $workdir/data

#TransDecoder.LongOrfs -t Trinity_cdhit90.fasta

#cd $workdir

## Step 1: BLAST processed Trinity contigs against Ceratopteris richardii, keep alignment info for analysis
##         and pull all unaligned ids to filter against trinity seqs file

# BLAST ORFs against Cr only - restrict query cover to >80% (at least 80% of trinity contig matches)
#diamond blastx --db ./ref/blastdb/nr.dmnd --query $workdir/data/Trinity_cdhit90.fasta.transdecoder_dir/longest_orfs.cds \
# --out ./data/7_blast/Trinity_cdhit90_500_blast_Cr_only.blastx \
# --outfmt 6 qseqid sseqid qlen slen evalue pident mismatch staxids stitle \
# --max-target-seqs 1 --evalue 1e-10 --threads 30 --taxonlist 49495 --unal 1 \
# --query-cover 80


#load required modules 
#module purge;
#module load bluebear
#module load bear-apps/2021b
#module load seqtk/1.3-GCC-11.2.0


# count how many trinity contigs did not match to Ceratopteris richardii 
#awk -F'\t' '$9 == "*" {count++} END {print "A total of " count " sequences did not have a significant match to the Ceratopteris richardii genome"}' ./data/7_blast/Trinity_cdhit90_500_blast_Cr_only.blastx > missing_genes_log.txt

# Pull any seqs that didn't match to Ceratopteris 
#awk -F'\t' '$9 == "*"' ./data/7_blast/Trinity_cdhit90_500_blast_Cr_only.blastx | cut -f 1 > ./data/7_blast/non_cr_500_filtered_ORFS.txt


# Filter trinity seqs by non_cr_500_filtered_ORFS
#seqtk subseq $workdir/data/Trinity_cdhit90.fasta.transdecoder_dir/longest_orfs.cds ./data/7_blast/non_cr_500_filtered_ORFS.txt > ./data/7_blast/Trinity_non_cr_500_ORFs.fasta

## Step 2: BLAST filtered Trinity contigs that did not match to Ceratopteris against the full nr database
# blast ORFs against nr database

#module purge;
#module load bluebear
#module load bear-apps/2022b
#module load DIAMOND/2.1.8-GCC-12.2.0

#diamond blastx --db ./ref/blastdb/nr.dmnd --query ./data/7_blast/Trinity_non_cr_500_ORFs.fasta \
# --out ./data/7_blast/Trinity_cdhit90_500_ORFS_nr_hits.blastx \
# --outfmt 6 qseqid sseqid qlen slen evalue pident mismatch staxids sscinames sskingdoms skingdoms sphylums stitle \
# --max-target-seqs 1 --evalue 1e-10 --threads 30

#load required modules 
#module purge;
#module load bluebear
#module load bear-apps/2021b
#module load seqtk/1.3-GCC-11.2.0

#Pull all sequences matching to "Viridiplantae", blast cr only and pull non hits -> missing genes 
#awk -F'\t' '$11 == "Viridiplantae" && $9 != "Ceratopteris richardii"' ./data/7_blast/Trinity_cdhit90_500_ORFS_nr_hits.blastx | cut -f 1 > ./data/7_blast/Trinity_cdhit90_500_ORFS_nr_hits_noncr_plant_ids.txt

# Filter Trinity_non_cr_500_ORFs by noncr_plant ids 
#seqtk subseq ./data/7_blast/Trinity_non_cr_500_ORFs.fasta ./data/7_blast/Trinity_cdhit90_500_ORFS_nr_hits_noncr_plant_ids.txt > ./data/7_blast/Trinity_cdhit90_500_ORFs_nocr_plant.fasta

## Step 3: Double check they dont match to Cr 
# BLAST sequences against Cr only without restricting query cover to >80%

#module purge;
#module load bluebear
#module load bear-apps/2022b
#module load DIAMOND/2.1.8-GCC-12.2.0

# BLAST ORFs against Cr only 
#diamond blastx --db ./ref/blastdb/nr.dmnd --query ./data/7_blast/Trinity_cdhit90_500_ORFs_nocr_plant.fasta \
# --out ./data/7_blast/Trinity_cdhit90_500_ORFs_plant_no_cr_cr.blastx \
# --outfmt 6 qseqid sseqid qlen slen evalue pident mismatch staxids stitle \
# --max-target-seqs 1 --evalue 1e-10 --threads 30 --taxonlist 49495 --unal 1 
# --query-cover 80


#load required modules 
module purge;
module load bluebear
module load bear-apps/2021b
module load seqtk/1.3-GCC-11.2.0

# Pull all seqs that did not match Ceratopteris -> 

# count how many trinity contigs did not match to Ceratopteris richardii 
#awk -F'\t' '$9 == "*" {count++} END {print "A total of " count " sequences did not have a significant match to the Ceratopteris richardii genome"}' ./data/7_blast/Trinity_cdhit90_500_ORFs_plant_no_cr_cr.blastx

# Pull any seqs that didn't match to Ceratopteris 
#awk -F'\t' '$9 == "*"' ./data/7_blast/Trinity_cdhit90_500_ORFs_plant_no_cr_cr.blastx | cut -f 1 > ./data/7_blast/missing_genes_ids.txt


# Filter trinity seqs by non_cr_500_filtered_ORFS
#seqtk subseq ./data/7_blast/Trinity_cdhit90_500_ORFs_nocr_plant.fasta ./data/7_blast/missing_genes_ids.txt > ./data/7_blast/missing_genes_salmon.fasta

# filter kallisto counts matrix 
awk 'NR==FNR {ids[$1]; next} FNR==1 || $1 in ids' ./data/7_blast/missing_genes_ids.txt ./data/6_rrna_filtered/kallisto.isoform.counts.matrix > ./data/7_blast/missing_genes_salmon_counts_matrix.txt

echo "kallisto gene counts matrix filtered by missing genes and saved to missing_gene_counts_matrix.txt"





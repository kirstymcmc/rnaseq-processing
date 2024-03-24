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

module purge;
module load bluebear
module load bear-apps/2022b
module load DIAMOND/2.1.8-GCC-12.2.0

# BLAST trinity contigs against nt database using DIAMOND
# To retrieve taxonomy information (needed to classify alignments by taxon), taxon info
# must be supplied in the form of a taxon map file. This file is a tab-delimited file
# These can be downloaded from NCBI's taxonomy database
# Files are in ref/blastdb/

workdir="/rds/projects/p/plackarg-spl-bioinformatic-analysis/full_dataset"
#cd $workdir/ref/blastdb

##make custom nr database with taxon ids 
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
#unzip taxdmp.zip
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz

#diamond makedb --in nr.gz -d nr --taxonmap prot.accession2taxid.gz --taxonnodes nodes.dmp --taxonnames names.dmp --threads 30


cd $workdir
mkdir -p ./data/7_blast

#filter trinity contigs to remove seqeunces less than 500bp in length
#perl ./rnaseq-processing/filter_fasta.pl ./data/Trinity_cdhit90.fasta 500 > ./data/Trinity_cdhit90_500.fasta


#diamond blastx --db ./ref/blastdb/nr.dmnd --query ./data/Trinity_cdhit90_500.fasta \
# --out ./data/7_blast/Trinity_cdhit90_500_nr.blastx \
# --outfmt 6 qseqid sseqid qstart qend sstart send evalue pident mismatch staxids sscinames sskingdoms skingdoms sphylums \
# --max-target-seqs 1 --evalue 1e-10 --threads 30


#load required modules 
module purge;
module load bluebear
module load bear-apps/2021b
module load seqtk/1.3-GCC-11.2.0

# Filter the Trinity_cdhit90_500_nr.blastx file for trinity contigs that are of the kingdom "Viridiplantae" but not
# the species Ceratopteris richardii 

awk -F'\t' '$13 == "Viridiplantae" && $11 != "Ceratopteris richardii"' ./data/7_blast/Trinity_cdhit90_500_nr.blastx | cut -f 1 > ./data/7_blast/non_cr_plant_500_filtered_contigs_ids.txt


# filter trinity seqs by no_cr_plant_filtered_contig_ids.txt and save to new file 

seqtk subseq ./data/Trinity_cdhit90_500.fasta ./data/7_blast/non_cr_plant_500_filtered_contigs_ids.txt > ./data/7_blast/Trinity_non_cr_plant_500_contigs.fasta


# run DIAMOND again but restrict to Ceratopteris richardii 
# purge and load modules
module purge;
module load bluebear
module load bear-apps/2022b
module load DIAMOND/2.1.8-GCC-12.2.0

diamond blastx --db ./ref/blastdb/nr.dmnd --query ./data/7_blast/Trinity_non_cr_plant_500_contigs.fasta \
 --out ./data/7_blast/Trinity_cdhit90_500_nr_nocr.blastx \
 --outfmt 6 qseqid sseqid qstart qend sstart send evalue pident mismatch staxids sscinames sskingdoms skingdoms sphylums \
 --taxonlist 49495 \
 --max-target-seqs 1 --evalue 1e-10 --threads 30 --unal 1 

# count how many trinity contigs did not match to Ceratopteris richardii 
awk -F'\t' '$11 == "*" {count++} END {print "A total of " count " sequences did have a significant match to the Ceratopteris richardii genome"}' ./data/7_blast/Trinity_cdhit90_500_nr_nocr.blastx

# pull all trinity "genes" that did not produce a hit to Ceratopteris richardii 

awk -F'\t' '$11 == "*"' ./data/7_blast/Trinity_cdhit90_500_nr_nocr.blastx | cut -f 1 > ./data/7_blast/missing_gene_ids.txt

echo "Missing gene IDs extracted and stored in missing_gene_ids.txt"

# filter trinity seqs by non-cr plant ids and save to file for downstream analysis 
#load required modules 
module purge;
module load bluebear
module load bear-apps/2021b
module load seqtk/1.3-GCC-11.2.0

seqtk subseq ./data/Trinity_cdhit90_500.fasta ./data/7_blast/missing_gene_ids.txt > ./data/7_blast/missing_genes.fasta

echo "Missing gene sequences extracted and stored in missing_genes.fasta"

# filter kallisto counts matrix 
awk 'NR==FNR {ids[$1]; next} FNR==1 || $1 in ids' ./data/7_blast/missing_gene_ids.txt ./data/6_rrna_filtered/kallisto.isoform.counts.matrix > ./data/7_blast/missing_genes_counts_matrix.txt

echo "kallisto gene counts matrix filtered by missing genes and saved to missing_gene_counts_matrix.txt"














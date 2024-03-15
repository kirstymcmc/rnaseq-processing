#!/bin/bash

#SBATCH --qos bbdefault
#SBATCH --account plackarg-spl-bioinformatic-analysis
#SBATCH --ntasks 30 # request 30 cores - necessary for parallel step
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 5-00:00 # this requests 5 days
#SBATCH --mail-type ALL 

# Reads from the preprocessing pipeline are in 6_rrna_filtered folder
# Trinity assembly is in $workdir folder
# Trinity tools are in script folder folder
# Define $pattern for automatic listing of samples
# trinity_sample_table.txt is sample file used with Trinity
# Requires dplyr, data.table and stringr (v1.5.0) packages installed in R

pattern="Np_"
workdir="../data"


################################################
# Removal of redundancy and quality estimation #
################################################
#Load required modules 
module purge;
module load bluebear
module load bear-apps/2022b
module load CD-HIT/4.8.1-GCC-12.2.0
module load Perl/5.36.0-GCCcore-12.2.0
module load Bowtie2/2.5.1-GCC-12.2.0

# Removing redundancy
echo "Removing redundancy from the assembly" | tee $workdir/postprocess_log.txt
cd-hit-est \
	-i $workdir/trinity_out_dir.Trinity.fasta \
	-o $workdir/Trinity_cdhit90.fasta -c 0.90 -n 9 -d 0 -M 0 -T 30 -s 0.9 -aS 0.9

echo "Redundancy removed, proceeding to basic statistics" | tee -a $workdir/postprocess_log.txt

# Basic statistics
#./TrinityStats.pl $workdir/Trinity_cdhit90.fasta > $workdir/Trinity_cdhit90_stats.txt
apptainer exec -e $workdir/6_rrna_filtered/trinityrnaseq.v2.15.1.simg /usr/local/bin/util/TrinityStats.pl trinity_out_dir.Trinity.fasta > Trinity_assembly.metrics

# Mapping reads to the assembly

## Creating bowtie2 db

bowtie2-build $workdir/Trinity_cdhit90.fasta \
	$workdir/Trinity.db --threads 30


## Mapping the reads
bowtie2 -p 12 -q --no-unal -k 20 \
-x '$workdir/Trinity.db' \
-1 '$workdir/trinity_out_dir/insilico_read_normalization/left.norm.fq' \
-2 '$workdir/trinity_out_dir/insilico_read_normalization/right.norm.fq' \
2> '$workdir/align_stats.txt' | samtools view -@10 -Sb -o bowtie2.bam

# Preparing new gene_trans_map for cdhit90.fasta
# Extracting sequence names
awk 'sub(/^>/, "")' $workdir/Trinity_cdhit90.fasta \
> Trinity_cdhit90_headers.txt

# Removing the unwanted part
awk '{print $1}' $workdir/Trinity_cdhit90_headers.txt \
> Trinity_cdhit90_header_filtered.txt

# Filtering the original gene_trans_map
# Use -v to inverse the selection
grep -Fwf $workdir/Trinity_cdhit90_header_filtered.txt \
$workdir/trinity_out_dir.Trinity.fasta.gene_trans_map \
> Trinity_cdhit90.fasta.gene_trans_map

# Kallisto
apptainer exec -e $workdir/6_rrna_filtered/trinityrnaseq.v2.15.1.simg /usr/local/bin/util/align_and_estimate_abundance.pl \
--transcripts $workdir/Trinity_cdhit90.fasta \
--gene_trans_map $workdir/Trinity_cdhit90.fasta.gene_trans_map \
--seqType fq --samples_file $workdir/trinity_sample_table.txt \
--est_method kallisto --aln_method bowtie2 --SS_lib_type RF \
--thread_count 24 --trinity_mode --prep_reference --output_dir $workdir/kallisto_outdir

# for Kallisto use this
find $pattern* -name "abundance.tsv" | tee quant_files.list

apptainer exec -e $workdir/6_rrna_filtered/trinityrnaseq.v2.15.1.simg /usr/local/bin/util/abundance_estimates_to_matrix.pl \
--est_method kallisto \
--gene_trans_map $workdir/Trinity_cdhit90.fasta.gene_trans_map \
--quant_files $workdir/quant_files.list --name_sample_by_basedir

# ExN50
apptainer exec -e $workdir/6_rrna_filtered/trinityrnaseq.v2.15.1.simg /usr/local/bin/util/misc/contig_ExN50_statistic.pl \
$workdir/kallisto.isoform.TMM.EXPR.matrix \
$workdir/Trinity_cdhit90.fasta transcript | tee ExN50.transcript.stats

apptainer exec -e $workdir/6_rrna_filtered/trinityrnaseq.v2.15.1.simg /usr/local/bin/util/misc/plot_ExN50_statistic.Rscript \
$workdir/ExN50.transcript.stats

xpdf $workdir/ExN50_plot.pdf

echo "Basic statistics finished. All output files generated." | tee -a $workdir/postprocess_log.txt
echo "Starting annotations using Trinotate" | tee -a $workdir/postprocess_log.txt
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
module load Biopython/1.81-foss-2022b

from Bio import SeqIO
from collections import defaultdict

# Function to parse GFF file and return annotations
def parse_gff(gff_file):
    annotations = defaultdict(lambda: defaultdict(list))
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
            annotations[seqid][feature_type].append((int(start), int(end), attributes))
    return annotations

# Function to infer introns from exon annotations
def infer_introns(annotations):
    introns = defaultdict(list)
    for seqid, features in annotations.items():
        exons = sorted(features["exon"], key=lambda x: x[0])  # Sort exons by start position
        for i in range(len(exons) - 1):
            intron_start = exons[i][1] + 1
            intron_end = exons[i + 1][0] - 1
            if intron_end > intron_start:  # Checking to ensure intron length is positive
                introns[seqid].append((intron_start, intron_end))
    return introns

# Function to infer intergenic regions from gene annotations
def infer_intergenic(annotations):
    intergenic_regions = defaultdict(list)
    for seqid, features in annotations.items():
        genes = sorted(features["gene"], key=lambda x: x[0])  # Sort genes by start position
        for i in range(len(genes) - 1):
            intergenic_start = genes[i][1] + 1
            intergenic_end = genes[i + 1][0] - 1
            if intergenic_end > intergenic_start:  # Checking to ensure region length is positive
                intergenic_regions[seqid].append((intergenic_start, intergenic_end))
    return intergenic_regions

# Function to save regions (introns or intergenic) to a BED file
#def save_regions_to_bed(regions, filename):
#    with open(filename, 'w') as f:
#        for seqid, region_list in regions.items():
#            for start, end in region_list:
                # Note: BED format is 0-based, but GFF is 1-based
                # Adjust start position by -1 for 0-based start position in BED
 #               f.write(f"{seqid}\t{start-1}\t{end}\n")
#

# Main
gff_file = "../ref/blastdb/Crichardii_676_v2.1.gene_exons.gff3"
annotations = parse_gff(gff_file)
introns = infer_introns(annotations)
intergenic_regions = infer_intergenic(annotations)

# Example: print first 5 inferred introns for Chr01
print("First 5 inferred introns for Chr01:", introns["Chr01"][:5])
# Example: print first 5 inferred intergenic regions for Chr01
print("First 5 inferred intergenic regions for Chr01:", intergenic_regions["Chr01"][:5])

# After calculating introns and intergenic regions
#introns_bed_file = "introns.bed"
#intergenic_regions_bed_file = "intergenic_regions.bed"
#save_regions_to_bed(introns, introns_bed_file)
#save_regions_to_bed(intergenic_regions, intergenic_regions_bed_file)
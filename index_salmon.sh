#!/bin/bash

#SBATCH --qos bbdefault
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 2:00:00 # this requests 2 hours (dd-hh:mm:ss)
#SBATCH --mem-per-cpu=6750M
#SBATCH --ntasks=54
#SBATCH --account plackarg-spl-bioinformatic-analysis
#SBATCH --mail-type ALL

module purge;
module load bluebear
module load bear-apps/2022b
module load Salmon/1.10.1-GCC-12.2.0



# This script creates an index for salmon using code from Salmon 1.10.2 documentation
# Requires both a reference genome and transcriptome in fasta format

# Download the reference transcriptome and genome for Ceratopteris richardii

# the latest genome annotation (Crichardii_676_v2.0.fa.gz) and transcriptome 
# (Crichardii_676_v2.1.cds.fa.gz) were downloaded from Phytozome 
# (https://data.jgi.doe.gov/refine-download/phytozome?organism=Crichardii&_gl=1*uimcht*_ga*MTQwODY4Mzk2NC4xNzA0OTgzMjI1*_ga_YBLMHYR3C2*MTcwOTAzNzAwNC42LjAuMTcwOTAzNzAwNC4wLjAuMA..&expanded=Phytozome-676)
# and moved to ./ref


# Extract names of genome targets 

#grep "^>" <(gunzip -c Crichardii_676_v2.0.fa.gz) | cut -d " " -f 1 > decoys.txt
#sed -i.bak -e 's/>//g' decoys.txt

#concatenate transcriptome and genome reference files 
#NOTE transcriptome targets must come first 

#cat Crichardii_676_v2.1.cds.fa.gz Crichardii_676_v2.0.fa.gz > gentrome.fa.gz

cd ../ref
#Now run salmon indexing 

salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index
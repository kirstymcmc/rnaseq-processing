#!/bin/bash

#SBATCH --qos bbdefault
#SBATCH --account plackarg-spl-bioinformatic-analysis
#SBATCH --mem-per-cpu=6750M
#SBATCH --ntasks 54 # request 54 cores for the job. 
#SBATCH --nodes 1 # restrict the job to a single node. Necessary if requesting more than --ntasks=1
#SBATCH --time 12:00:00 # this requests 12 hours, 
#SBATCH --mail-type ALL

set -e

module purge;
module load bluebear
module load bear-apps/2022b
module load Salmon/1.10.1-GCC-12.2.0



# Define work folders 
workdir="../data"
fqdir="$workdir/2_trimmed"

refdir=../ref/salmon_index
outdir="$workdir/3_quants"

echo "Salmon quantification started" >> $workdir/pipeline_log.txt
echo Current Date and Time is: `date +"%Y-%m-%d %T"` >> $workdir/pipeline_log.txt

# Check if the output directory exists, if not, create it
if [ ! -d "$outdir" ]; then
  mkdir -p "$outdir"
fi

for r1 in "$fqdir"/*1.fq.gz; do
    r2=${r1%1.fq.gz}2.fq.gz
    out_name=$(basename "$r1" | cut -d '_' -f 1-3)
    if [[ -f $r2 ]] ; then

salmon quant -i $refdir -l A \
-1 "$r1" -2 "$r2" \
-p 54 --validateMappings -o $outdir/$out_name --writeUnmappedNames \
--gcBias --seqBias
    else
        echo "$r2 not found" >&2
    fi
done


echo "Salmon quantification completed. Extracting unmapped reads..." >> $workdir/pipeline_log.txt
echo Current Date and Time is: `date +"%Y-%m-%d %T"` >> $workdir/pipeline_log.txt
# Unmapped names are written to 3_quants/${f}/aux_info/unmapped_names.txt
# extract names and seqs using extract_unmapped.sh
# takes as input (in order): 1. data folder (directory containing 1_raw, 2_trimmed etc)
#                            2. fastq folder containing all reads (eg 2_trimmed)
#                            3. output folder name (created within the script
#                            4. flag (one of "u", "d", "m1", "m2") 

./extract_unmapped.sh $workdir $fqdir $outdir "u"

echo "Unmapped reads extracted" >> $workdir/pipeline_log.txt
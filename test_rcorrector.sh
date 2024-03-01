module purge;
module load bluebear
module load bear-apps/2021b
module load Python/3.9.6-GCCcore-11.2.0
module load Rcorrector/1.0.5-GCC-11.2.0

# Removing erroneous k-mers from Illumina paired-end reads
# Run with the highest possible number of cores
SECONDS=0
echo "Running Rcorrector on raw reads" | tee -a $workdir/pipeline_log.txt
mkdir 5_corrected
echo "Corrected reads will be in 5_corrected folder" | tee -a $workdir/pipeline_log.txt

# Raw reads are in folder '4_unmapped'
fqdir=$workdir/4_unmapped
for r1 in "$fqdir"/*1_unmapped.fastq; do
    r2=${r1%1_unmapped.fastq.gz}2_unmapped.fastq
    if [[ -f $r2 ]] ; then
        run_rcorrector.pl -1 $r1 -2 $r2 -t 24 -od '$workdir/5_corrected/'
    else
        echo "$r2 not found" >&2
    fi
done
echo "Time needed to finish Rcorrector step on raw reads: $SECONDS seconds" >> $workdir/pipeline_log.txt
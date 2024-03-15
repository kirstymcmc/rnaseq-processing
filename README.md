
# RNAseq processing with Salmon, extraction of unmapped reads and subsequent de novo transcriptome synthesis 

This collection of scripts takes raw fastq files from RNAseq, aligns them to the transcriptome 
using salmon, extracts read pairs that did not map, preprocesses and assembles them de novo using the 
Trinity pipeline. 

Steps of the workflow are as follows: 

## Preprocess raw reads
Quality of raw reads is assessed before and after trimming with script salmon_run.sh
### Initial quality control
FastQC and MultiQC are used to assess quality of raw reads 
### Trimming adapters and low quality bases 
[Trim Galore](https://github.com/FelixKrueger/TrimGalore) is used for adapter and quality trimming. 
### QC trimmed reads
FastQC and multiQC are used to assess read quality after trimming. 

## Quantify reads with [Salmon](https://github.com/COMBINE-lab/salmon)
### Generate Salmon index
Use index_salmon.sh to generate an index for Ceratopteris richardii
### Quantify reads and extract unmapped
Run Salmon quantification of trimmed reads with quantify_salmon.sh. After quantification, unmapped reads (where both reads in a pair did not map) are extracted and written to 4_unmapped.

## Preprocess unmapped reads 
All preprocessing steps are in the file preprocess_unmapped.sh.
Code roughly follows pipeline from (https://github.com/matevzl533/Noccaea_praecox_transcriptome/tree/main)
### Initial quality control
FastQC and MultiQC are used to assess quality of "raw" unmapped reads
### Removing erroneous k-mers from Illimina paired-end reads 
[rCorrector](https://github.com/mourisl/Rcorrector) is used to tag reads in the fastq output as corrected or uncorrectable. rcorrector is a tool specifically designed for kmer-bases read error correction of RNA-seq data.
### Discard read pairs for which one or both reads is deemed unfixable
Uses a python script from the Harvard Informatics GitHub repository [TranscriptomeAssemblyTools](https://github.com/harvardinformatics/TranscriptomeAssemblyTools). The script has been updated to Python3.
### Remove unwanted rrna reads with Bowtie2
From Silva, the SSUParc and LSUParc fasta files were downloaded (https://ftp.arb-silva.de/?pk_vid=8352a8ccf0ead1d7168388545541b6c1). Before running bowtie2-build, SSUParc and LSUParc were concatenated and U translated to T. 
```shell
cat *.fasta > SILVA.db
awk '/^[^>]/ { gsub(/U/,"T"); print; next }1' SILVA.db > SLVA.db
```
### Run fastqc on processed reads 
Re-run QC from step 1.

## de novo assemble with Trinity
### Make sample table text file
Trinity accepts a text file via --samples_file rather than looping through reads [see here](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity). Run make_sample_table.py and provide the directory containing your clean reads. 
### Run trinity
[Trinity](https://github.com/trinityrnaseq/trinityrnaseq) is used for de novo transcriptome assembly with default parameters. Script to run Trinity is in trinity.sh.   

## Post-processing 
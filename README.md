RNAseq processing with Salmon, extraction of unmapped reads and subsequent de novo 
transcriptome synthesis. 

This collection of scripts takes raw fastq files from RNAseq, aligns them to the transcriptome 
using salmon, extracts read pairs that did not map and de novo assembles them using the 
Trinity pipeline. 

Steps of the workflow are as follows: 

•add QC steps 
1. Generate salmon index with index_salmon.sh
2. quantify reads with quantify_salmon.sh
3. 

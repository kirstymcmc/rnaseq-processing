#!/bin/bash

# Check usage
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 data_folder fastq_folder output_folder flag"
    exit 1
fi

data_folder="$1"
fastq_folder="$2"
output_folder="$3"
flag="$4"  # The flag to filter by, e.g., "u"

# Create a directory for the unmapped names files if it doesn't exist
mkdir -p "$data_folder/unmapped_names"

# Check if the unmapped_names folder is already populated
num_files=$(find "$data_folder/unmapped_names" -type f | wc -l)
if [ "$num_files" -eq 0 ]; then
    echo "Moving and renaming unmapped_names.txt files..."
    FILES="$data_folder/3_quants/*"
    for f in $FILES; do 
        f_name=$(basename "$f")
        f_name=${f_name%_quant}
        mv "$f/aux_info/unmapped_names.txt" "$data_folder/unmapped_names/${f_name}.txt"
    done
else
    echo "The unmapped_names folder is already populated. Skipping moving and renaming step."
fi

# Ensure the output folder for the FASTQ files exists
mkdir -p "$output_folder"

# Process each unmapped_names file in the unmapped folder
for unmapped_file in "$data_folder/unmapped_names"/*.txt; do
    base_name=$(basename "$unmapped_file" .txt)

    # Filter to keep only lines with the specified flag
    grep " $flag" "$unmapped_file" | cut -d ' ' -f 1 > "${unmapped_file%.txt}_filtered.txt"
    filtered_unmapped="${unmapped_file%.txt}_filtered.txt"

    # Loop through both pairs of FASTQ files
    for i in 1 2; do
        input_file="$fastq_folder/${base_name}_${i}_trimmo.fq.gz"
        output_file="$output_folder/${base_name}_${i}_unmapped.fastq"

        # Process FASTQ file
        zcat "$input_file" | paste - - - - | awk -v unmapped="$filtered_unmapped" 'BEGIN{while((getline < unmapped) > 0) ids["@"$1];}{if($1 in ids) print $0;}' | tr '\t' '\n' > "$output_file"
    done
done

import os
import json

def parse_meta_info(json_file_path):
    """Parse the meta_info.json file to extract required information."""
    with open(json_file_path, 'r') as file:
        data = json.load(file)
    # Extracting required information
    return {
        "num_processed": data.get("num_processed", 0),
        "num_mapped": data.get("num_mapped", 0),
        "num_decoy_fragments": data.get("num_decoy_fragments", 0),
    }

def parse_unmapped_names(unmapped_names_file_path):
    """Parse the unmapped names file to count occurrences of each flag."""
    counts = {"u": 0, "d": 0, "m1": 0, "m2": 0}
    with open(unmapped_names_file_path, 'r') as file:
        for line in file:
            flag = line.strip().split()[-1]
            if flag in counts:
                counts[flag] += 1
    return counts

def summarize_samples(main_directory_path):
    # Header for the output table
    print("Sample\tTotal Reads\tMapped\tDecoys\tUnmapped\tDovetail (Discordant)\tm1\tm2")

    for root, dirs, files in os.walk(main_directory_path):
        for sample_dir in dirs:
            json_path = os.path.join(root, sample_dir, "aux_info", "meta_info.json")
            unmapped_names_path = os.path.join(root, sample_dir, "logs", "unmapped_names.txt")  # Adjust if the path differs

            if os.path.exists(json_path):
                meta_info = parse_meta_info(json_path)
                
            if os.path.exists(unmapped_names_path):
                unmapped_counts = parse_unmapped_names(unmapped_names_path)

            # Adjust the print statement according to your structure; this is just an example
            print(f"{sample_dir}\t{meta_info['num_processed']}\t{meta_info['num_mapped']}\t{meta_info['num_decoy_fragments']}\t{unmapped_counts['u']}\t{unmapped_counts['d']}\t{unmapped_counts['m1']}\t{unmapped_counts['m2']}")

main_directory_path = '/path/to/your/main/directory'  # Update this path
summarize_samples(main_directory_path)

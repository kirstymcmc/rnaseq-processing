
from Bio import SeqIO
from collections import defaultdict
import re

def parse_gff(gff_file):
    annotations = defaultdict(lambda: defaultdict(list))
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            seqid, _, feature_type, start, end, _, _, _, attributes = parts
            annotations[seqid][feature_type].append((int(start), int(end), attributes))
    return annotations


def parse_gff_for_gene_exon_hierarchy(gff_file_path):
    """
    Parses a GFF file and returns a hierarchy of genes, their transcripts, and associated exons.

    Parameters:
    gff_file_path (str): Path to the GFF file.

    Returns:
    dict: A nested dictionary with gene IDs as keys, each containing a dictionary of transcript IDs,
          each of which contains a list of exon (start, end) tuples.
    """
    genes = defaultdict(lambda: defaultdict(list))

    with open(gff_file_path, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue  # Skip header lines
            parts = line.strip().split("\t")
            feature_type = parts[2]
            attributes = parts[8]

            # More robust attribute parsing that accounts for complex scenarios
            attr_dict = {match.group(1): match.group(2) for match in re.finditer(r'(\w+)=([^;]+)', attributes)}

            if feature_type == "exon":
                parent_id = attr_dict.get("Parent")
                # Use regex to extract the full gene ID from the parent transcript ID
                gene_id_match = re.match(r'(.*?\.\d+)', parent_id)
                if gene_id_match:
                    gene_id = gene_id_match.group(1)
                else:
                    gene_id = "Unknown"  # or some default handling if no match

                # Now store the exon start and end positions, instead of just the ID
                start, end = int(parts[3]), int(parts[4])
                genes[gene_id][parent_id].append((start, end))

    return genes



def infer_introns_within_genes(gene_exon_hierarchy):
    introns = defaultdict(list)
    for gene_id, transcripts in gene_exon_hierarchy.items():
        for transcript_id, exons in transcripts.items():
            # Ensure exons are sorted by their start positions
            sorted_exons = sorted(exons, key=lambda x: x[0])  # Sort exons by start position
            for i in range(len(sorted_exons) - 1):
                # Calculate intron start and end based on exon positions
                intron_start = sorted_exons[i][1] + 1  # End of current exon + 1
                intron_end = sorted_exons[i + 1][0] - 1  # Start of next exon - 1
                if intron_end > intron_start:  # Ensure intron length is positive
                    introns[gene_id].append((intron_start, intron_end))
    return introns


def infer_intergenic(annotations):
    intergenic_regions = defaultdict(list)
    for seqid, features in annotations.items():
        genes = sorted(features["gene"], key=lambda x: x[0])
        for i in range(len(genes) - 1):
            intergenic_start = genes[i][1] + 1
            intergenic_end = genes[i + 1][0] - 1
            if intergenic_end > intergenic_start:
                intergenic_regions[seqid].append((intergenic_start, intergenic_end))
    return intergenic_regions

def output_to_bed(features_data, output_file_path, feature_label):
    """
    Write genomic features data to a BED file.

    Parameters:
    - features_data: A dictionary where keys are sequence IDs and values are lists of tuples, each tuple representing
      the start and end positions of a feature.
    - output_file_path: Path to the output BED file.
    - feature_label: A label indicating the type of feature (e.g., "exon", "intron", "intergenic") to include in the BED file.
    """
    with open(output_file_path, 'w') as out_file:
        for seqid, features in features_data.items():
            for start, end in features:
                # Adjust start position by -1 for BED format, which is 0-based
                bed_start = start - 1
                bed_end = end
                # BED format: chrom, chromStart, chromEnd, name, score, strand
                # Assuming '+' strand for all features; adjust accordingly if strand information is available
                out_file.write(f"{seqid}\t{bed_start}\t{bed_end}\t{seqid}_{feature_label}\t0\t+\n")


def main():
    gff_file_path = "../../ref/blastdb/Crichardii_676_v2.1.gene_exons.gff3"
    
    annotations = parse_gff(gff_file_path)
    gene_exon_hierarchy = parse_gff_for_gene_exon_hierarchy(gff_file_path)
    intron_data = infer_introns_within_genes(gene_exon_hierarchy)
    intergenic_data = infer_intergenic(annotations)

    # Prepare exon data for all seqids
    exon_data = {}
    for seqid in annotations:
        if "exon" in annotations[seqid]:
            exon_data_for_seqid = [(start, end) for start, end, _ in annotations[seqid]["exon"]]
            if exon_data_for_seqid:  # If there are exons for this seqid
                exon_data[seqid] = exon_data_for_seqid

    # Output to single files
    output_to_bed(exon_data, 'exons.bed', 'exon')
    output_to_bed(intron_data, 'introns.bed', 'intron')
    output_to_bed(intergenic_data, 'intergenic.bed', 'intergenic')

    print("Consolidated BED files for exons, introns, and intergenic regions have been generated.")


if __name__ == "__main__":
    main()

# Main
#gff_file = "../../ref/blastdb/Crichardii_676_v2.1.gene_exons.gff3"
#annotations = parse_gff(gff_file)
#gene_exon_hierarchy = parse_gff_for_gene_exon_hierarchy(gff_file)
#introns_within_genes = infer_introns_within_genes(gene_exon_hierarchy)
#intergenic_regions = infer_intergenic(annotations)

# Debugging Outputs

#print("First 5 inferred intergenic regions for Chr01:", intergenic_regions.get("Chr01", [])[:5])
#introns_within_genes = infer_introns_within_genes(gene_exon_hierarchy)

# Debug: Print first few inferred introns for a specific gene
#for gene_id, intron_positions in list(introns_within_genes.items())[:1]:  # Adjust as needed
#    print(f"Gene ID: {gene_id}, First few introns: {intron_positions[:5]}")


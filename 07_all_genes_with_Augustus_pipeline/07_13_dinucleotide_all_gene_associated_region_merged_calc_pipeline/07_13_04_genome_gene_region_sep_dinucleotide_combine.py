import pandas as pd

# Define file paths
file1_path = "/home/y-kondo/protostome_lifespan_model/Genome_dinucleotide_count_Github_check/genome_dinucleotide_summary_with_mono_rates_Github_check.txt"

file2_path = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/dinucleotide_count_merged/summary_dinucleotide_rates_gene_associated_-2kb_+1kb_renamed_Github_check.csv"

output_path = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/dinucleotide_count_merged/combined_dinucleotide_genome_gene_associated_-2kb_+1kb_Github_check.csv"

# Read the files
file1 = pd.read_csv(file1_path, sep="\t")
file2 = pd.read_csv(file2_path, sep=",")

# Add "Region" column to file1 and set it to "genome"
file1["Region"] = "genome"

# Rename columns in file1
file1.rename(columns={
    "Genome_Name": "Species"
}, inplace=True)


# Get the unique species in each file
species_file1 = set(file1["Species"])
species_file2 = set(file2["Species"])

# Find species in file2 but not in file1
species_only_in_file2 = species_file2 - species_file1

# Output the results
print(f"Species in file2 but not in file1: {len(species_only_in_file2)}")
print(species_only_in_file2)

# Find species in file1 but not in file2
species_only_in_file1 = species_file1 - species_file2

# Output the results
print(f"Species in file1 but not in file2: {len(species_only_in_file1)}")
print(species_only_in_file1)

# Keep only species common to both files
common_species = species_file1 & species_file2
file1_filtered = file1[file1["Species"].isin(common_species)].copy()
file2_filtered = file2[file2["Species"].isin(common_species)].copy()

# Combine the filtered data
combined_data = pd.concat([file1_filtered, file2_filtered], ignore_index=True)

# Save the combined file
combined_data.to_csv(output_path, index=False)

print(f"Combined file (only species present in both files) saved to {output_path}")

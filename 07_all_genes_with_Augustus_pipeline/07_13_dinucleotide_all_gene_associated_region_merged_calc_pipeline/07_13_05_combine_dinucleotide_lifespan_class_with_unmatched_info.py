import pandas as pd

# File paths
dinucleotide_rates_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/dinucleotide_count_merged/combined_dinucleotide_genome_gene_associated_-2kb_+1kb_Github_check.csv"
lifespan_file = "/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/busco_species_lifespan_with_class_20250507.csv"

output_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/dinucleotide_count_merged/combined_dinucleotide_lifespan_class_Github_check.csv"
unmatched_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/dinucleotide_count_merged/unmatched_species_Github_check.csv"

# Load the data
dinucleotide_rates_df = pd.read_csv(dinucleotide_rates_file)
lifespan_df = pd.read_csv(lifespan_file)

# Merge the data on 'Species' and 'Assembly Identifier'
combined_df = pd.merge(
    dinucleotide_rates_df,
    lifespan_df,
    left_on="Species",
    right_on="Assembly Identifier",
    how="left"
)

# Extract unmatched species names
unmatched_species = combined_df.loc[combined_df["Organism Name"].isna(), "Species"].dropna().drop_duplicates()

# Save only the list of unmatched species names
unmatched_species.to_csv(unmatched_file, index=False, header=["Unmatched_Species"])
print(f"{len(unmatched_species)} unmatched species saved to {unmatched_file}")

# Remove unmatched rows (where Organism Name is NaN)
combined_df = combined_df[~combined_df["Organism Name"].isna()].copy()

# Drop the 'Species' column since it's the same as 'Assembly Identifier'
combined_df = combined_df.drop(columns=["Species"])

# Save the combined data
combined_df.to_csv(output_file, index=False)
print(f"Combined data saved to {output_file}")
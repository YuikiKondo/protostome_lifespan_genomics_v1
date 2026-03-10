import pandas as pd

# Input files
counts_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/gene_associated_-2kb_+1kb_merged_total_bp_per_species.csv"

main_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/dinucleotide_count_merged/dinucleotide_PGLS/filtered_species_with_zero_lifecycle_6_categories_Github_check_with_intergenic_counts_rates_OE_with_tip_name.csv"

# Output file
out_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/dinucleotide_count_merged/dinucleotide_PGLS/filtered_species_with_zero_lifecycle_6_categories_Github_check_with_intergenic_counts_rates_OE_with_gene_associated_total_bp.csv"

# Load
counts_df = pd.read_csv(counts_file)
main_df = pd.read_csv(main_file)

# Keep only what we need (and ensure Species is unique)
counts_df = counts_df[["Species", "Total_bp_gene_associated_region_including_N"]].drop_duplicates(subset=["Species"])

# Merge (left join keeps all rows in main_df)
merged = main_df.merge(
    counts_df,
    how="left",
    left_on="Assembly Identifier",
    right_on="Species"
)

# Drop the extra Species column created by the merge
merged = merged.drop(columns=["Species"])

# Safety check
n_missing = merged["Total_bp_gene_associated_region_including_N"].isna().sum()
if n_missing > 0:
    print(f"WARNING: {n_missing} rows did not get a value (no match).")

# Save
merged.to_csv(out_file, index=False)
print("Wrote:", out_file)

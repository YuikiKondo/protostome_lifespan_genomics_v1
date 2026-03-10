import pandas as pd

# Define file paths
mechanics_file = "/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/Downstream_analysis_Github_check/all_species_sliding_window_downstream_dinucleotide_counts_stepsize_10.csv"
match_file = "/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/TSS_Github_check/Assembly_species_match.csv"    # Use the same file as TSS sliding window analysis
output_file = "/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/Downstream_analysis_Github_check/all_species_sliding_window_downstream_dinucleotide_counts_stepsize_10_with_species.csv"

# Load datasets
mechanics_df = pd.read_csv(mechanics_file)
match_df = pd.read_csv(match_file)

# Merge dataframes based on "Assembly" column
merged_df = mechanics_df.merge(match_df, on="Assembly", how="left")

# Save the updated dataframe
merged_df.to_csv(output_file, index=False)

print(f"Updated file saved as: {output_file}")

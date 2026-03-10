import pandas as pd

# Input and output files
input_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/up_downstreams_range_refined_exon_intron_separated/dinucleotide_count/combined_dinucleotide_lifespan_class_Github_check.csv"
output_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/up_downstreams_range_refined_exon_intron_separated/dinucleotide_count/filtered_dinucleotide_lifespan_class_Github_check.csv"

## Exclude 2 species with no BUSCOs hit and 1 species with only 1 hit gene
## Exclude 5 species for which AUGUSTUS failed to run as well
exclude_species_file = "/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/tetranucleotide_count/exclude_species_list_busco_zero_hit_and_augustus_failed_species.csv"

# Read the exclusion list from the CSV file
exclude_species_df = pd.read_csv(exclude_species_file)
species_to_exclude = exclude_species_df.iloc[:, 0].tolist()  # Assuming species names are in the first column

# Read the main data file
data = pd.read_csv(input_file)

# Exclude specified species
filtered_data = data[~data["Organism Name"].isin(species_to_exclude)]

# Save the filtered data
filtered_data.to_csv(output_file, index=False)

print(f"Filtered data saved to {output_file}")

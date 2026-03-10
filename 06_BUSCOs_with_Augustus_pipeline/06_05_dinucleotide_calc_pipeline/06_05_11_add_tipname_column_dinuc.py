import pandas as pd

# Define file paths 
busco_file = "/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_exon_intron_separated/gene_region_separated_sequences_no_internal_stop_3divisible/tetranucleotide_count/DNA_mechanics/PGLS/busco_species_tree_tip_alternative_names.csv"

filtered_species_file = "/home/y-kondo/protostome_lifespan_model_augustus/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/dinucleotide_count/filtered_species_with_zero_lifecycle_6_categories_Github_check.csv"

output_file = "/home/y-kondo/protostome_lifespan_model_augustus/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/dinucleotide_count/dinucleotide_PGLS/filtered_species_with_tip_name_dinuc_Github_check.csv"

# Load the CSV files
busco_df = pd.read_csv(busco_file)
filtered_species_df = pd.read_csv(filtered_species_file)

# Rename 'Gene_Length_Including_N' to 'Genome_size_bp'
if 'Gene_Length_Including_N' in filtered_species_df.columns:
    filtered_species_df = filtered_species_df.rename(columns={'Gene_Length_Including_N': 'Genome_size_bp'})
else:
    print("Warning: 'Gene_Length_Including_N' column not found in the filtered_species file.")

# Normalize blanks to NaN (in case some blanks are empty strings)
filtered_species_df['Genome_size_bp'] = (
    filtered_species_df['Genome_size_bp']
    .replace(r'^\s*$', pd.NA, regex=True)
)

# Build a species → genome size mapping from the 'genome' region rows
if {'Organism Name', 'Region', 'Genome_size_bp'}.issubset(filtered_species_df.columns):
    genome_rows = filtered_species_df.loc[
        filtered_species_df['Region'].astype(str).str.lower() == 'genome',
        ['Organism Name', 'Genome_size_bp']
    ].dropna(subset=['Genome_size_bp'])

    # Optional: check for conflicting genome sizes within the same species
    conflicts = (
        genome_rows.groupby('Organism Name')['Genome_size_bp']
        .nunique(dropna=True)
        .reset_index()
    )
    conflicts = conflicts[conflicts['Genome_size_bp'] > 1]
    if not conflicts.empty:
        print("Warning: Conflicting 'Genome_size_bp' values found for these species:")
        print(conflicts['Organism Name'].tolist())

    # Prefer the first non-null value per species
    size_map = (
        genome_rows.drop_duplicates(subset=['Organism Name'], keep='first')
        .set_index('Organism Name')['Genome_size_bp']
    )

    # Fill missing genome sizes across all regions for the same species
    filtered_species_df['Genome_size_bp'] = filtered_species_df['Genome_size_bp'].fillna(
        filtered_species_df['Organism Name'].map(size_map)
    )
else:
    print("Warning: Missing required columns among {'Organism Name','Region','Genome_size_bp'}.")

# Merge based on "Organism Name" (filtered_species_df) and "Species" (busco_df)
merged_df = filtered_species_df.merge(
    busco_df[['Species', 'Tip_Name']],
    left_on="Organism Name",
    right_on="Species",
    how="left"
)

# Drop the redundant "Species" column
if 'Species' in merged_df.columns:
    merged_df = merged_df.drop(columns=['Species'])

# Save the merged file
merged_df.to_csv(output_file, index=False)
print(f"File saved successfully: {output_file}")

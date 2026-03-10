import pandas as pd

# Load the input CSV file
input_file = "filtered_species_with_zero_lifecycle_6_categories_Github_check.csv"
df = pd.read_csv(input_file)

# Specify columns to keep
columns_to_keep = ['Organism Name', 'Assembly Identifier', 'Region', 'Sum_of_valid_dinucleotide_counts']
columns_to_keep += [col for col in df.columns if ('_Rate' in col or '_Count' in col) and col not in ['A_Rate', 'T_Rate', 'G_Rate', 'C_Rate', 'A_Count', 'T_Count', 'G_Count', 'C_Count']]

# Filter the dataframe
df_filtered = df[columns_to_keep]

# Order rows by Region and then Organism Name
df_filtered = df_filtered.sort_values(by=['Organism Name', 'Region'])

# Save the filtered data to a new CSV file with your specified name
output_file = "di_rates_and_counts_with_lifespan_Github_check.csv"
df_filtered.to_csv(output_file, index=False)

print(f"Filtered data saved to {output_file}")

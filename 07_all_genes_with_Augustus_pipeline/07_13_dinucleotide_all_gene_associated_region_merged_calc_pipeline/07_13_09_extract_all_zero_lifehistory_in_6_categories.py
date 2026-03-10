import pandas as pd

# Input and output file paths
input_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/dinucleotide_count_merged/filtered_dinucleotide_lifespan_class_Github_check.csv"
output_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/dinucleotide_count_merged/filtered_species_with_zero_lifecycle_6_categories_Github_check.csv"

# Load the dataset
data = pd.read_csv(input_file)

# Filter the dataset to include only species where all 5 categories are 0
filtered_data = data[
    (data["Parthenogenesis"] == 0) &
    (data["Fission Budding"] == 0) &
    (data["Highly dormant eggs"] == 0) &
    (data["Dormant state in other than egg stage"] == 0) &
    (data["More than one type of life cycles"] == 0) &
    (data["Eusociality"] == 0)
]

# Save the filtered dataset to a CSV file
filtered_data.to_csv(output_file, index=False)

print(f"Filtered data saved to: {output_file}")

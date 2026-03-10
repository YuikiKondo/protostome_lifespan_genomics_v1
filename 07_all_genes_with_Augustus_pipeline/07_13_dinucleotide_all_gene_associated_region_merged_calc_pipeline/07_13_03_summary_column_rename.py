import pandas as pd

# Define input and output file paths
input_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/dinucleotide_count_merged/summary_dinucleotide_rates_gene_associated_-2kb_+1kb.csv"

output_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/dinucleotide_count_merged/summary_dinucleotide_rates_gene_associated_-2kb_+1kb_renamed_Github_check.csv"

# Read the CSV file
df = pd.read_csv(input_file)

# Rename columns
df.rename(columns={"Total_Valid_Dinucleotide_Counts": "Sum_of_valid_dinucleotide_counts"}, inplace=True)
df.columns = [col.replace(" Count", "_Count").replace(" Rate", "_Rate").replace(" O/E", "_OE") for col in df.columns]

# Save the modified CSV file
df.to_csv(output_file, index=False)

print("Column names updated and saved to", output_file)
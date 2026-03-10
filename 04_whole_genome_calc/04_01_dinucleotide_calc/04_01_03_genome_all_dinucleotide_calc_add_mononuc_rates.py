import pandas as pd

# Define input and output file paths
input_file = "/home/y-kondo/protostome_lifespan_model/Genome_dinucleotide_count_Github_check/genome_dinucleotide_summary_all_Github_check.txt"
output_file = "/home/y-kondo/protostome_lifespan_model/Genome_dinucleotide_count_Github_check/genome_dinucleotide_summary_with_mono_rates_Github_check.txt"

# Read the input file
df = pd.read_csv(input_file, sep="\t")

# Calculate total nucleotide count
df["Total_Nucleotides"] = df["A_Count"] + df["T_Count"] + df["G_Count"] + df["C_Count"]

# Calculate rates
df["A_Rate"] = df["A_Count"] / df["Total_Nucleotides"]
df["T_Rate"] = df["T_Count"] / df["Total_Nucleotides"]
df["G_Rate"] = df["G_Count"] / df["Total_Nucleotides"]
df["C_Rate"] = df["C_Count"] / df["Total_Nucleotides"]

# Drop the temporary total nucleotide column
df.drop(columns=["Total_Nucleotides"], inplace=True)

# Save the updated file
df.to_csv(output_file, sep="\t", index=False)

print(f"Updated file saved as: {output_file}")
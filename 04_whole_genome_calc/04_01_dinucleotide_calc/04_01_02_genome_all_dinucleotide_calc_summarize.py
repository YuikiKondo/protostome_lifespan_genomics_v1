import os

# Directory containing the files
directory = "/home/y-kondo/protostome_lifespan_model/Genome_dinucleotide_count_Github_check"

# List of files to concatenate
files = [
    "genome_dinucleotide_summary0_100.txt",
    "genome_dinucleotide_summary100_200.txt",
    "genome_dinucleotide_summary200_300.txt",
    "genome_dinucleotide_summary300_400.txt",
    "genome_dinucleotide_summary400_500.txt",
    "genome_dinucleotide_summary500_650.txt",
]

# Output file name
output_file = "genome_dinucleotide_summary_all_Github_check.txt"

# Full path to the output file
output_path = os.path.join(directory, output_file)

# Concatenating files
with open(output_path, "w") as outfile:
    for idx, fname in enumerate(files):
        file_path = os.path.join(directory, fname)
        with open(file_path, "r") as infile:
            if idx == 0:
                # Write the header from the first file
                outfile.write(infile.read())
            else:
                # Skip header for all subsequent files
                next(infile)  # Skip the first line
                outfile.write(infile.read())

print(f"Concatenation completed: {output_path}")

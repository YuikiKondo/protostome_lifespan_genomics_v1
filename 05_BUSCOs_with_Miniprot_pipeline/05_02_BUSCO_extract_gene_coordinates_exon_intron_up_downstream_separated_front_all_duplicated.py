import os
import pandas as pd

# Define base directory containing species directories
base_busco_dir = '/home/y-kondo/protostome_lifespan_model/BUSCO_directory'
# Define output directory for all processed files
output_base_dir = '/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/busco_coordinates_including_dup'
os.makedirs(output_base_dir, exist_ok=True)  # Ensure the output directory exists

# Directory names to search for within each species folder
busco_subdirs = ['single_copy_busco_sequences', 'multi_copy_busco_sequences', 'fragmented_busco_sequences']

# Function to process .gff files in a specified directory and output to a file
def process_gff_files(input_dir, output_filename):
    all_coordinates = []

    # Loop through each .gff file in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".gff"):
            gene_id = filename.split(".")[0]  # Extract gene ID from filename
            gff_file_path = os.path.join(input_dir, filename)

            mRNAs = {}
            exons_by_mRNA = {}

            # Read the file and gather mRNA and exon data
            with open(gff_file_path, 'r') as file:
                for line in file:
                    if line.startswith("#"):
                        continue

                    fields = line.strip().split('\t')
                    if len(fields) < 9:
                        continue  # skip incomplete lines

                    reference = fields[0]
                    feature = fields[2]
                    start = int(fields[3]) - 1
                    end = int(fields[4])
                    score = float(fields[5]) if fields[5] != "." else 0
                    strand = fields[6]
                    attributes = fields[8]

                    if feature == 'mRNA':
                        mRNA_id = attributes.split("ID=")[1].split(";")[0]
                        mRNAs[mRNA_id] = (score, strand)

                    elif feature == 'CDS':
                        parent_id = attributes.split("Parent=")[1].split(";")[0]
                        if parent_id not in exons_by_mRNA:
                            exons_by_mRNA[parent_id] = []
                        exons_by_mRNA[parent_id].append((reference, start, end, score, strand))

            # Process all mRNAs instead of just the best one with the highest hit score
            for mRNA_id in mRNAs:
                score, strand = mRNAs[mRNA_id]
                exons = exons_by_mRNA.get(mRNA_id, [])
                if not exons:
                    continue

                exons = sorted(exons, key=lambda x: x[1], reverse=(strand == '-'))

                TSS_start = exons[0][1] if strand == '+' else exons[0][2]

                # Upstream regions
                if strand == '+':
                    upstream1_start = TSS_start - 300
                    upstream1_end = TSS_start - 100
                    upstream2_start = TSS_start - 100
                    upstream2_end = TSS_start + 100
                else:
                    upstream1_start = TSS_start + 100
                    upstream1_end = TSS_start + 300
                    upstream2_start = TSS_start - 100
                    upstream2_end = TSS_start + 100

                all_coordinates.append([f"{gene_id}_{mRNA_id}_upstream1", exons[0][0], upstream1_start, upstream1_end, strand, score])
                all_coordinates.append([f"{gene_id}_{mRNA_id}_upstream2", exons[0][0], upstream2_start, upstream2_end, strand, score])

                # Exons and inferred introns
                for i, (reference, start, end, exon_score, strand) in enumerate(exons):
                    exon_number = f"{gene_id}_{mRNA_id}_{i+1}_exon"
                    all_coordinates.append([exon_number, reference, start, end, strand, exon_score])

                    if i < len(exons) - 1:
                        if strand == '+':
                            intron_start = end
                            intron_end = exons[i + 1][1]
                        else:
                            intron_start = exons[i + 1][2]
                            intron_end = start
                        if intron_start <= intron_end:
                            intron_number = f"{gene_id}_{mRNA_id}_{i+1}_intron"
                            all_coordinates.append([intron_number, reference, intron_start, intron_end, strand, exon_score])

                # Downstream regions
                if strand == '+':
                    downstream1_start = exons[-1][2] - 100
                    downstream1_end = exons[-1][2]
                    downstream2_start = exons[-1][2]
                    downstream2_end = exons[-1][2] + 200
                else:
                    downstream1_start = exons[-1][1]
                    downstream1_end = exons[-1][1] + 100
                    downstream2_start = exons[-1][1] - 200
                    downstream2_end = exons[-1][1]

                all_coordinates.append([f"{gene_id}_{mRNA_id}_downstream1", exons[-1][0], downstream1_start, downstream1_end, strand, exons[-1][3]])
                all_coordinates.append([f"{gene_id}_{mRNA_id}_downstream2", exons[-1][0], downstream2_start, downstream2_end, strand, exons[-1][3]])

    columns = ["Region", "Reference", "Start", "End", "Strand", "Score"]
    coordinates_df = pd.DataFrame(all_coordinates, columns=columns)
    coordinates_df.to_csv(output_filename, sep='\t', index=False)
    print(f"Data saved to: {output_filename}")

# Loop through each species directory
for species_dir in os.listdir(base_busco_dir):
    species_path = os.path.join(base_busco_dir, species_dir)
    run_dir = os.path.join(species_path, 'run_metazoa_odb10', 'busco_sequences')
    if os.path.isdir(run_dir):
        for subdir in busco_subdirs:
            input_path = os.path.join(run_dir, subdir)
            if os.path.isdir(input_path):
                output_file = os.path.join(
                    output_base_dir,
                    f"{species_dir}_{subdir.replace('single_copy_busco_sequences', 'complete_genes').replace('multi_copy_busco_sequences', 'duplicated_genes').replace('fragmented_busco_sequences', 'fragmented_genes')}.txt"
                )
                process_gff_files(input_path, output_file)

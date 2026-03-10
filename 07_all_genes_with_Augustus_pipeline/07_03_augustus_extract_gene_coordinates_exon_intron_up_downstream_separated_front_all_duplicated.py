import os
import re
import pandas as pd

# Define base directories
augustus_gff_dir = '/home/y-kondo/protostome_lifespan_model_augustus/augustus_predictions/completed_augustus_final_dir'
output_base_dir = '/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/up_downstreams_range_refined_exon_intron_separated/augustus_coordinates'
os.makedirs(output_base_dir, exist_ok=True)

# Utility to extract attribute values
def extract_attribute(attributes_str, key):
    match = re.search(rf'{key} "([^"]+)"', attributes_str)
    return match.group(1) if match else None

# Process a single .gff file in Augustus format (modified for whole genome gff)
def process_augustus_gff_whole_genome(gff_file_path):
    all_coordinates = []
    transcripts = {}

    with open(gff_file_path, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            reference = fields[0]
            feature = fields[2]
            start = int(fields[3]) - 1
            end = int(fields[4])
            score = float(fields[5]) if fields[5] != "." else 0
            strand = fields[6]
            attributes = fields[8]

            if feature == 'transcript':
                transcript_id = attributes.strip() if ";" not in attributes else extract_attribute(attributes, 'transcript_id')
                gene_id = extract_attribute(attributes, 'gene_id') or transcript_id
                transcripts[transcript_id] = {
                    'gene_id': gene_id,
                    'reference': reference,
                    'strand': strand,
                    'score': score,
                    'exons': []
                }

            elif feature == 'CDS':
                transcript_id = extract_attribute(attributes, 'transcript_id')
                if transcript_id in transcripts:
                    transcripts[transcript_id]['exons'].append((reference, start, end, score, strand))

    for transcript_id, info in transcripts.items():
        gene_id = info['gene_id']
        strand = info['strand']
        score = info['score']
        exons = sorted(info['exons'], key=lambda x: x[1], reverse=(strand == '-'))

        if not exons:
            continue

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

        all_coordinates.append([f"{gene_id}_{transcript_id}_upstream1", info['reference'], upstream1_start, upstream1_end, strand, score])
        all_coordinates.append([f"{gene_id}_{transcript_id}_upstream2", info['reference'], upstream2_start, upstream2_end, strand, score])

        for i, (ref, start, end, exon_score, strand) in enumerate(exons):
            exon_name = f"{gene_id}_{transcript_id}_{i+1}_exon"
            all_coordinates.append([exon_name, ref, start, end, strand, exon_score])

            if i < len(exons) - 1:
                if strand == '+':
                    intron_start = end
                    intron_end = exons[i + 1][1]
                else:
                    intron_start = exons[i + 1][2]
                    intron_end = start
                if intron_start <= intron_end:
                    intron_name = f"{gene_id}_{transcript_id}_{i+1}_intron"
                    all_coordinates.append([intron_name, ref, intron_start, intron_end, strand, exon_score])

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

        all_coordinates.append([f"{gene_id}_{transcript_id}_downstream1", exons[-1][0], downstream1_start, downstream1_end, strand, exons[-1][3]])
        all_coordinates.append([f"{gene_id}_{transcript_id}_downstream2", exons[-1][0], downstream2_start, downstream2_end, strand, exons[-1][3]])

    return all_coordinates

# Main loop over each .gff file
for gff_filename in os.listdir(augustus_gff_dir):
    if not gff_filename.endswith('.gff'):
        continue

    species_id = gff_filename.replace('_augustus.gff', '')
    gff_path = os.path.join(augustus_gff_dir, gff_filename)

    coords = process_augustus_gff_whole_genome(gff_path)

    output_file = os.path.join(output_base_dir, f"{species_id}_augustus_all_genes.txt")
    df = pd.DataFrame(coords, columns=["Region", "Reference", "Start", "End", "Strand", "Score"])
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Data saved to: {output_file} ({'with data' if not df.empty else 'empty'})")

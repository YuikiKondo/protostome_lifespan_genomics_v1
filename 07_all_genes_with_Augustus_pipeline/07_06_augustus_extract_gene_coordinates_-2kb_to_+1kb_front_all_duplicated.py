import os
import re
import pandas as pd

# Define base directories
augustus_gff_dir = '/home/y-kondo/protostome_lifespan_model_augustus/augustus_predictions/completed_augustus_final_dir'
output_base_dir = '/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_coordinates'
os.makedirs(output_base_dir, exist_ok=True)

# Utility to extract attribute values from Augustus GFF (key "value")
def extract_attribute(attributes_str, key):
    m = re.search(rf'{re.escape(key)}\s+"([^"]+)"', attributes_str)
    return m.group(1) if m else None

def process_augustus_gff_whole_genome(gff_file_path, upstream_bp=2000, downstream_bp=1000, clamp_start=True):
    """
    For each transcript, compute a single region:
      gene_associated_-2kb_+1kb = [ATG-2000, STOP+1000] (strand-aware)
    using CDS coordinates.
    Output coordinates are 0-based, end-exclusive (consistent with bedtools).
    """
    transcripts = {}  # transcript_id -> dict

    with open(gff_file_path, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split('\t')
            if len(fields) < 9:
                continue

            reference = fields[0]
            feature = fields[2]
            start = int(fields[3]) - 1  # 0-based
            end = int(fields[4])        # end-exclusive
            score = float(fields[5]) if fields[5] != "." else 0.0
            strand = fields[6]
            attributes = fields[8]

            if feature == 'transcript':
                # Augustus sometimes stores transcript_id as an attribute; sometimes the full attributes string is an ID-like blob
                transcript_id = extract_attribute(attributes, 'transcript_id')
                if transcript_id is None:
                    transcript_id = attributes.strip()

                gene_id = extract_attribute(attributes, 'gene_id') or transcript_id

                transcripts.setdefault(transcript_id, {
                    'gene_id': gene_id,
                    'reference': reference,
                    'strand': strand,
                    'score': score,
                    'cds': []
                })

            elif feature == 'CDS':
                transcript_id = extract_attribute(attributes, 'transcript_id')
                if transcript_id is None:
                    # If transcript_id is missing, we can’t safely attach this CDS to a transcript
                    continue

                if transcript_id not in transcripts:
                    # Create a placeholder transcript entry if CDS appears before transcript line
                    gene_id = extract_attribute(attributes, 'gene_id') or transcript_id
                    transcripts[transcript_id] = {
                        'gene_id': gene_id,
                        'reference': reference,
                        'strand': strand,
                        'score': score,
                        'cds': []
                    }

                transcripts[transcript_id]['cds'].append((reference, start, end, score, strand))

    all_coordinates = []

    for transcript_id, info in transcripts.items():
        cds_list = info['cds']
        if not cds_list:
            continue

        strand = info['strand']
        gene_id = info['gene_id']
        ref_set = {x[0] for x in cds_list}

        # If a transcript spans multiple contigs (rare), skip to avoid wrong extraction
        if len(ref_set) != 1:
            continue

        reference = next(iter(ref_set))

        cds_min_start = min(x[1] for x in cds_list)
        cds_max_end = max(x[2] for x in cds_list)

        # Start/stop codon positions on genome (strand-aware)
        if strand == '+':
            start_codon_pos = cds_min_start
            stop_codon_pos = cds_max_end
            window_start = start_codon_pos - upstream_bp
            window_end = stop_codon_pos + downstream_bp
        else:
            start_codon_pos = cds_max_end
            stop_codon_pos = cds_min_start
            window_start = stop_codon_pos - downstream_bp
            window_end = start_codon_pos + upstream_bp

        # Clamp
        if clamp_start:
            window_start = max(0, window_start)

        # Basic sanity
        if window_end <= window_start:
            continue

        region_name = f"{gene_id}_{transcript_id}_gene_associated_-2kb_+1kb"
        all_coordinates.append([region_name, reference, window_start, window_end, strand, info['score']])

    return all_coordinates

# Main loop over each .gff file
for gff_filename in os.listdir(augustus_gff_dir):
    if not gff_filename.endswith('.gff'):
        continue

    species_id = gff_filename.replace('_augustus.gff', '').replace('.gff', '')
    gff_path = os.path.join(augustus_gff_dir, gff_filename)

    coords = process_augustus_gff_whole_genome(gff_path, upstream_bp=2000, downstream_bp=1000, clamp_start=True)

    output_file = os.path.join(output_base_dir, f"{species_id}_augustus_all_genes_gene_associated_-2kb_+1kb.txt")
    df = pd.DataFrame(coords, columns=["Region", "Reference", "Start", "End", "Strand", "Score"])
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Data saved to: {output_file} ({'with data' if not df.empty else 'empty'})")

## Purpouse: extract gene coordinates with each gene region separately (upstream1, upstream2, exon, intron, downstream1, downstream2) using gff files outputted by running BUSCO with Augustus. The coordinates of Complete/Duplicated/Fragmented genes are each outputted in a separate file by referring to busco_gene_status_summary.csv (such as /home/username/protostome_lifespan_model_augustus/BUSCO_directory/BUSCO_all_species/BUSCO_GCA_029721355.1_KUTeg_genomic.fna/run_metazoa_odb10/busco_sequences/busco_gene_status_summary.csv).

## example output:
#(base) [username@usernumber busco_coordinates_including_dup]$ ls|head -n 6
#BUSCO_GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna_complete_genes.txt
#BUSCO_GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna_duplicated_genes.txt
#BUSCO_GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna_fragmented_genes.txt
#BUSCO_GCA_000002075.2_AplCal3.0_genomic.fna_complete_genes.txt
#BUSCO_GCA_000002075.2_AplCal3.0_genomic.fna_duplicated_genes.txt
#BUSCO_GCA_000002075.2_AplCal3.0_genomic.fna_fragmented_genes.txt
#(base) [username@usernumber busco_coordinates_including_dup]$ head #BUSCO_GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna_complete_genes.txt
#Region	Reference	Start	End	Strand	Score
#102804at33208_r2.m1.g2.t1_upstream1	AE014296.5	1876278	1876478	-	0.08
#102804at33208_r2.m1.g2.t1_upstream2	AE014296.5	1876078	1876278	-	0.08
#102804at33208_r2.m1.g2.t1_1_exon	AE014296.5	1876109	1876178	-	0.82
#102804at33208_r2.m1.g2.t1_1_intron	AE014296.5	1875915	1876109	-	0.82
#102804at33208_r2.m1.g2.t1_2_exon	AE014296.5	1875651	1875915	-	1.0
#102804at33208_r2.m1.g2.t1_2_intron	AE014296.5	1875586	1875651	-	1.0
#102804at33208_r2.m1.g2.t1_3_exon	AE014296.5	1875286	1875586	-	0.92
#102804at33208_r2.m1.g2.t1_3_intron	AE014296.5	1875218	1875286	-	0.92
#102804at33208_r2.m1.g2.t1_4_exon	AE014296.5	1874040	1875218	-	1.0

import os
import re
import pandas as pd

base_busco_dir = '/home/username/protostome_lifespan_model_augustus/BUSCO_directory/BUSCO_all_species'
output_base_dir = '/home/username/protostome_lifespan_model_augustus/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/busco_coordinates_including_dup'
os.makedirs(output_base_dir, exist_ok=True)

# Utility to extract attribute values like transcript_id or gene_id
def extract_attribute(attributes_str, key):
    match = re.search(rf'{key}\s+"([^"]+)"', attributes_str)
    return match.group(1) if match else None

# Parse the gene status summary file and return a dictionary
def parse_gene_status(status_file):
    df = pd.read_csv(status_file)
    gene_status_map = {}
    for _, row in df.iterrows():
        if row['Single_copy'] == 1:
            gene_status_map[row['GeneID']] = 'complete_genes'
        elif row['Multi_copy'] == 1:
            gene_status_map[row['GeneID']] = 'duplicated_genes'
        elif row['Fragmented'] == 1:
            gene_status_map[row['GeneID']] = 'fragmented_genes'
    return gene_status_map

# Process a single .gff file in Augustus format
def process_augustus_gff(gff_file_path, gene_id, transcript_id):
    all_coordinates = []

    transcripts = {}
    exons = []

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
                transcripts[transcript_id] = (score, strand)

            elif feature == 'CDS':
                exons.append((reference, start, end, score, strand))

    if transcript_id not in transcripts or not exons:
        return all_coordinates  # skip if no valid transcript/exon

    score, strand = transcripts[transcript_id]
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

    all_coordinates.append([f"{gene_id}_{transcript_id}_upstream1", exons[0][0], upstream1_start, upstream1_end, strand, score])
    all_coordinates.append([f"{gene_id}_{transcript_id}_upstream2", exons[0][0], upstream2_start, upstream2_end, strand, score])

    for i, (reference, start, end, exon_score, strand) in enumerate(exons):
        exon_name = f"{gene_id}_{transcript_id}_{i+1}_exon"
        all_coordinates.append([exon_name, reference, start, end, strand, exon_score])

        if i < len(exons) - 1:
            if strand == '+':
                intron_start = end
                intron_end = exons[i + 1][1]
            else:
                intron_start = exons[i + 1][2]
                intron_end = start
            if intron_start <= intron_end:
                intron_name = f"{gene_id}_{transcript_id}_{i+1}_intron"
                all_coordinates.append([intron_name, reference, intron_start, intron_end, strand, exon_score])

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

# Main processing loop
for species_dir in os.listdir(base_busco_dir):
    species_path = os.path.join(base_busco_dir, species_dir)
    run_dir = os.path.join(species_path, 'run_metazoa_odb10')
    augustus_dir = os.path.join(run_dir, 'augustus_output/gff')
    status_file = os.path.join(run_dir, 'busco_sequences', 'busco_gene_status_summary.csv')

    if os.path.isdir(augustus_dir) and os.path.isfile(status_file):
        gene_status_map = parse_gene_status(status_file)

        status_output = {
            'complete_genes': [],
            'duplicated_genes': [],
            'fragmented_genes': []
        }

        for gff_filename in os.listdir(augustus_dir):
            if not gff_filename.endswith('.gff'):
                continue

            gene_id_raw = gff_filename.split('.')[0]
            if gene_id_raw not in gene_status_map:
                continue

            status = gene_status_map[gene_id_raw]
            gff_path = os.path.join(augustus_dir, gff_filename)

            transcript_id = None
            gene_id = gene_id_raw  # fallback

            with open(gff_path) as f:
                for line in f:
                    if "\ttranscript\t" in line:
                        fields = line.strip().split('\t')
                        attributes = fields[8]
                        # Case 1: plain transcript ID (no attributes)
                        if ';' not in attributes and '=' not in attributes and '\t' not in attributes:
                            transcript_id = attributes.strip()
                        else:
                            transcript_id = extract_attribute(attributes, 'transcript_id')
                            gene_id_extracted = extract_attribute(attributes, 'gene_id')
                            if gene_id_extracted:
                                gene_id = gene_id_extracted
                        break
                else:
                    continue  # skip if no transcript_id found

            coords = process_augustus_gff(gff_path, gene_id, transcript_id)
            status_output[status].extend(coords)

        # Output per gene status
        for status, records in status_output.items():
            output_file = os.path.join(output_base_dir, f"{species_dir}_{status}.txt")
            df = pd.DataFrame(records, columns=["Region", "Reference", "Start", "End", "Strand", "Score"])
            df.to_csv(output_file, sep='\t', index=False)
            print(f"{status} data saved to: {output_file} ({'with data' if not df.empty else 'empty'})")

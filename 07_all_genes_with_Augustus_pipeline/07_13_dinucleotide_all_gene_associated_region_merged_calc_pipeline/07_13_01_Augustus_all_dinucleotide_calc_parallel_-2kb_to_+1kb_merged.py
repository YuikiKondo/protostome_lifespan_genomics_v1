import os
import traceback
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
from collections import Counter

def calculate_dinucleotide_metrics(sequence):
    sequence = sequence.upper()
    length = len(sequence)

    bases = ['A', 'T', 'G', 'C']
    valid_dinucleotides = [a + b for a in bases for b in bases]

    dinucleotide_counts = Counter(
        sequence[i:i + 2] for i in range(length - 1)
        if sequence[i:i + 2] in valid_dinucleotides
    )
    valid_dinucleotide_counts = sum(dinucleotide_counts.values())

    nucleotide_counts = Counter(sequence)
    a_count = nucleotide_counts['A']
    t_count = nucleotide_counts['T']
    g_count = nucleotide_counts['G']
    c_count = nucleotide_counts['C']

    total_mono_count = a_count + t_count + g_count + c_count

    mono_rates = {
        base: (nucleotide_counts[base] / total_mono_count) if total_mono_count > 0 else 0
        for base in ['A', 'T', 'G', 'C']
    }

    expected_counts = {
        dn: (mono_rates[dn[0]] * mono_rates[dn[1]]) * valid_dinucleotide_counts
        if valid_dinucleotide_counts > 0 else 0
        for dn in valid_dinucleotides
    }

    oe_ratios = {
        dn: (dinucleotide_counts[dn] / expected_counts[dn]) if expected_counts[dn] > 0 else 0
        for dn in expected_counts
    }

    rates = {
        dn: (dinucleotide_counts[dn] / valid_dinucleotide_counts) if valid_dinucleotide_counts > 0 else 0
        for dn in valid_dinucleotides
    }

    return a_count, t_count, g_count, c_count, valid_dinucleotide_counts, dinucleotide_counts, oe_ratios, rates

def process_fasta_file(fasta_file, output_file):
    try:
        with open(output_file, "a") as out_f:
            for record in SeqIO.parse(fasta_file, "fasta"):
                header = record.description.split()

                # Expected Augustus header:
                # >geneid augustus scaffold start end strand score
                if len(header) < 7:
                    continue

                gene_id = header[0]
                source = header[1]          # "augustus"
                scaffold = header[2]
                gene_start = int(header[3])
                gene_end = int(header[4])
                strand = header[5]
                score = header[6]

                sequence = str(record.seq)
                gene_length_including_n = len(sequence)

                a_count, t_count, g_count, c_count, valid_dinucleotide_counts, dinucleotide_counts, oe_ratios, rates = \
                    calculate_dinucleotide_metrics(sequence)

                out_f.write(
                    f"{gene_id}\t{source}\t{scaffold}\t{score}\t{gene_start}\t{gene_end}\t{strand}\t"
                    f"{gene_length_including_n}\t{valid_dinucleotide_counts}\t{a_count}\t{t_count}\t{g_count}\t{c_count}"
                )

                for dn in [a + b for a in 'ATGC' for b in 'ATGC']:
                    count = dinucleotide_counts.get(dn, 0)
                    oe = oe_ratios.get(dn, 0)
                    rate = rates.get(dn, 0)
                    out_f.write(f"\t{count}\t{oe:.6f}\t{rate:.6f}")

                out_f.write("\n")

    except Exception as e:
        print(f"Error processing {fasta_file}: {e}")
        traceback.print_exc()

def process_species_files(species, species_files, output_file):
    try:
        with open(output_file, "w") as out_f:
            header = [
                "GeneID", "Source", "Scaffold", "Alignment_score", "Gene_Start", "Gene_End", "Strand",
                "Gene_Length_Including_N", "Valid_Dinucleotide_Counts", "A_Count", "T_Count", "G_Count", "C_Count"
            ]
            for dn in [a + b for a in 'ATGC' for b in 'ATGC']:
                header.extend([f"{dn}_Count", f"{dn}_OE", f"{dn}_Rate"])
            out_f.write("\t".join(header) + "\n")

        for fasta_file in species_files:
            process_fasta_file(fasta_file, output_file)

    except Exception as e:
        print(f"Error processing species {species}: {e}")
        traceback.print_exc()

def process_all_fasta_files(input_dir, output_dir, num_processes):
    os.makedirs(output_dir, exist_ok=True)

    suffix = "_augustus_all_genes_gene_associated_-2kb_+1kb_merged.fasta"

    species_files_map = {}
    for fasta_file in os.listdir(input_dir):
        if fasta_file.endswith(".fasta"):
            species = fasta_file
            if species.endswith(suffix):
                species = species[:-len(suffix)]
            else:
                # fallback: strip extension
                species = os.path.splitext(species)[0]

            fasta_path = os.path.join(input_dir, fasta_file)
            species_files_map.setdefault(species, []).append(fasta_path)

    print(f"Processing {len(species_files_map)} species with {num_processes} processes...")

    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = []
        for species, species_files in species_files_map.items():
            output_file = os.path.join(output_dir, f"{species}_gene_summary.txt")
            futures.append(executor.submit(process_species_files, species, species_files, output_file))

        for future in futures:
            future.result()

    print("Finished processing all species.")

if __name__ == "__main__":
    input_dir = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/augustus_sequences_merged/"
    output_dir = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/dinucleotide_count_merged/augustus_gene_summaries_Github_check/"

    num_processes = max(1, os.cpu_count() - 1)
    process_all_fasta_files(input_dir, output_dir, num_processes)

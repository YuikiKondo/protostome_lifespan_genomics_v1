import os
import sys
import traceback
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
from collections import Counter
import multiprocessing

def calculate_dinucleotide_metrics(sequence):
    """Calculate nucleotide counts, dinucleotide counts, and O/E ratios."""
    sequence = sequence.upper()
    length = len(sequence)

    # Generate all possible dinucleotides
    bases = ['A', 'T', 'G', 'C']
    valid_dinucleotides = [a + b for a in bases for b in bases]
    
    # Count overlapping dinucleotides
    dinucleotide_counts = Counter(sequence[i:i + 2] for i in range(length - 1) if sequence[i:i + 2] in valid_dinucleotides)
    valid_dinucleotide_counts = sum(dinucleotide_counts.values())

    # Count individual nucleotides
    nucleotide_counts = Counter(sequence)
    a_count = nucleotide_counts['A']
    t_count = nucleotide_counts['T']
    g_count = nucleotide_counts['G']
    c_count = nucleotide_counts['C']
    
    # Total mononucleotide count
    total_mono_count = sum(nucleotide_counts[base] for base in ['A', 'T', 'G', 'C'])

    # Mononucleotide frequencies
    mono_rates = {
        base: nucleotide_counts[base] / total_mono_count if total_mono_count > 0 else 0
        for base in ['A', 'T', 'G', 'C']
    }

    # Expected dinucleotide counts using mononucleotide frequencies
    expected_counts = {
        dn: (
            mono_rates[dn[0]] *
            mono_rates[dn[1]]
        ) * valid_dinucleotide_counts if valid_dinucleotide_counts > 0 else 0
        for dn in valid_dinucleotides
    }


    # Calculate O/E ratios
    oe_ratios = {dn: (dinucleotide_counts[dn] / expected_counts[dn]) if expected_counts[dn] > 0 else 0
                 for dn in expected_counts}

    # Calculate rates
    rates = {dn: (dinucleotide_counts[dn] / valid_dinucleotide_counts) if valid_dinucleotide_counts > 0 else 0
             for dn in valid_dinucleotides}

    return a_count, t_count, g_count, c_count, valid_dinucleotide_counts, dinucleotide_counts, oe_ratios, rates

def process_fasta_file(fasta_file, output_file, gene_type):
    """Process a single FASTA file and write gene metrics to the output file."""
    try:
        with open(output_file, "a") as out_f:
            for record in SeqIO.parse(fasta_file, "fasta"):
                header = record.description.split()
                busco_id = header[0]
                status = header[1]

                if gene_type == "missing":
                    out_f.write(f"{busco_id}\t{status}\n")
                else:
                    if len(header) < 6:
                        continue  # Skip incomplete headers
                    scaffold = header[2]
                    gene_start, gene_end = int(header[3]), int(header[4])
                    strand, score = header[5], header[6]
                    sequence = str(record.seq)

                    gene_length_including_n = len(sequence)
                    a_count, t_count, g_count, c_count, valid_dinucleotide_counts, dinucleotide_counts, oe_ratios, rates = calculate_dinucleotide_metrics(sequence)

                    out_f.write(f"{busco_id}\t{status}\t{scaffold}\t{score}\t{gene_start}\t{gene_end}\t{strand}\t{gene_length_including_n}\t{valid_dinucleotide_counts}\t{a_count}\t{t_count}\t{g_count}\t{c_count}")

                    for dinucleotide in [a + b for a in 'ATGC' for b in 'ATGC']:
                        count = dinucleotide_counts.get(dinucleotide, 0)
                        oe_ratio = oe_ratios.get(dinucleotide, 0)
                        rate = rates.get(dinucleotide, 0)
                        out_f.write(f"\t{count}\t{oe_ratio:.6f}\t{rate:.6f}")

                    out_f.write("\n")
    except Exception as e:
        print(f"Error processing {fasta_file}: {e}")
        traceback.print_exc()

def process_species_files(species, species_files, output_file):
    """Process all files for a given species in parallel."""
    try:
        with open(output_file, "w") as out_f:
            header = ["Busco ID", "Status", "Scaffold", "Alignment score", "Gene Start", "Gene End", "Strand", 
                      "Gene Length Including N", "Valid Dinucleotide Counts", "A Count", "T Count", 
                      "G Count", "C Count"]
            for dinucleotide in [a + b for a in 'ATGC' for b in 'ATGC']:
                header.extend([f"{dinucleotide} Count", f"{dinucleotide} O/E", f"{dinucleotide} Rate"])
            out_f.write("\t".join(header) + "\n")

        for fasta_file, gene_type in species_files:
            process_fasta_file(fasta_file, output_file, gene_type)

    except Exception as e:
        print(f"Error processing species {species}: {e}")
        traceback.print_exc()

def process_all_fasta_files(input_dir, output_dir, num_processes):
    """Process all FASTA files for all species using multiprocessing."""
    os.makedirs(output_dir, exist_ok=True)

    species_files_map = {}
    for fasta_file in os.listdir(input_dir):
        if fasta_file.endswith(".fasta"):
            species = fasta_file.split("_genomic.fna")[0]
            gene_type = fasta_file.split("_")[-2]  
            fasta_path = os.path.join(input_dir, fasta_file)
            species_files_map.setdefault(species, []).append((fasta_path, gene_type))

    print(f"Processing {len(species_files_map)} species with {num_processes} processes...")

    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = []
        for species, species_files in species_files_map.items():
            output_file = os.path.join(output_dir, f"{species}_gene_summary.txt")
            futures.append(executor.submit(process_species_files, species, species_files, output_file))

        for future in futures:
            future.result()  # Wait for completion

    print("Finished processing all species.")

if __name__ == "__main__":
    input_dir = "/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/busco_sequences/"
    output_dir = "/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/dinucleotide_count/busco_gene_summaries_Github_check/"

    num_processes = max(1, os.cpu_count() - 1)  # Uses all but one CPU core
    process_all_fasta_files(input_dir, output_dir, num_processes)

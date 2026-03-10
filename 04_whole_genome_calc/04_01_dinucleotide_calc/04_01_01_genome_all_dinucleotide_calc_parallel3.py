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
    valide_dinucleotide_counts = sum(dinucleotide_counts.values())

    # Count individual nucleotides
    nucleotide_counts = Counter(sequence)
    a_count = nucleotide_counts['A']
    t_count = nucleotide_counts['T']
    g_count = nucleotide_counts['G']
    c_count = nucleotide_counts['C']
    
    return a_count, t_count, g_count, c_count, valide_dinucleotide_counts, dinucleotide_counts

def write_header(output_file):
    """Write the header row if the file doesn't exist."""
    if not os.path.exists(output_file):  # Only write the header if the file is new
        dinucleotides = [a + b for a in 'ATGC' for b in 'ATGC']
        header = ["Genome_Name", "Gene_Length_Including_N", "Sum_of_valid_dinucleotide_counts", "A_Count", "T_Count", "G_Count", "C_Count"]
        for dn in dinucleotides:
            header.extend([f"{dn}_Count", f"{dn}_OE", f"{dn}_Rate"])
        
        with open(output_file, "w") as f:
            f.write("\t".join(header) + "\n")  # Tab-separated header row

def process_fasta_file(fasta_file, output_file):
    """Process a single genome FASTA file and write a single summary row for the entire genome."""
    try:
        genome_name = os.path.basename(fasta_file).replace("_genomic.fna", "")

        # Initialize cumulative values
        total_length_including_n = 0
        total_valide_dinucleotide_counts = 0
        total_nucleotide_counts = Counter()
        total_dinucleotide_counts = Counter()

        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence = str(record.seq)
            total_length_including_n += len(sequence)

            # Get nucleotide and dinucleotide counts for the sequence
            a_count, t_count, g_count, c_count, valide_dinucleotide_counts, dinucleotide_counts = calculate_dinucleotide_metrics(sequence)

            # Accumulate across sequences
            total_valide_dinucleotide_counts += valide_dinucleotide_counts
            total_nucleotide_counts.update({'A': a_count, 'T': t_count, 'G': g_count, 'C': c_count})
            total_dinucleotide_counts.update(dinucleotide_counts)

        valid_dinucleotides = [a + b for a in 'ATGC' for b in 'ATGC']

        # Total mono count
        total_mono_count = sum(total_nucleotide_counts[base] for base in ['A', 'T', 'G', 'C'])

        # Mononucleotide frequencies
        mono_rates = {
            base: total_nucleotide_counts[base] / total_mono_count if total_mono_count > 0 else 0
            for base in ['A', 'T', 'G', 'C']
        }

        # Expected dinucleotide counts using mono_rates
        expected_counts = {
            dn: (
                mono_rates[dn[0]] *
                mono_rates[dn[1]]
            ) * total_valide_dinucleotide_counts if total_valide_dinucleotide_counts > 0 else 0
            for dn in valid_dinucleotides
        }

        # O/E and Rate
        oe_ratios = {
            dn: (total_dinucleotide_counts[dn] / expected_counts[dn]) if expected_counts[dn] > 0 else 0
            for dn in valid_dinucleotides
        }
        rates = {
            dn: total_dinucleotide_counts[dn] / total_valide_dinucleotide_counts if total_valide_dinucleotide_counts > 0 else 0
            for dn in valid_dinucleotides
        }

        # Write the aggregated results as a single row
        with open(output_file, "a") as out_f:
            out_f.write(f"{genome_name}\t{total_length_including_n}\t{total_valide_dinucleotide_counts}\t"
                        f"{total_nucleotide_counts['A']}\t{total_nucleotide_counts['T']}\t{total_nucleotide_counts['G']}\t{total_nucleotide_counts['C']}")

            for dinucleotide in valid_dinucleotides:
                count = total_dinucleotide_counts.get(dinucleotide, 0)
                oe_ratio = oe_ratios.get(dinucleotide, 0)
                rate = rates.get(dinucleotide, 0)
                out_f.write(f"\t{count}\t{oe_ratio:.6f}\t{rate:.6f}")

            out_f.write("\n")

    except Exception as e:
        print(f"Error processing {fasta_file}: {e}")
        traceback.print_exc()


def process_all_fasta_files(input_dir, output_dir, num_processes, max_files=100):
    """Process up to max_files genome FASTA files using multiprocessing."""
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "genome_dinucleotide_summary200_300.txt")

    # Ensure header row is written before processing
    write_header(output_file)

    fasta_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith(".fna")][200:300]
    
    print(f"Processing {len(fasta_files)} genome files with {num_processes} processes...")

    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = []
        for fasta_file in fasta_files:
            futures.append(executor.submit(process_fasta_file, fasta_file, output_file))

        for future in futures:
            future.result()  # Wait for completion

    print("Finished processing selected genome files.")

if __name__ == "__main__":
    input_dir = "/home/y-kondo/protostome_lifespan_model/ProtostomeGenomes_unzipped/"
    output_dir = "/home/y-kondo/protostome_lifespan_model/Genome_dinucleotide_count_Github_check/"
    
    num_processes = max(1, os.cpu_count() - 1)  # Uses all but one CPU core
    process_all_fasta_files(input_dir, output_dir, num_processes, max_files=100)

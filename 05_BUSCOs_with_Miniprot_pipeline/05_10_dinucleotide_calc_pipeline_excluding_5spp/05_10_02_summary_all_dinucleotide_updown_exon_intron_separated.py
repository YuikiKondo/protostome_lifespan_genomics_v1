import os
import pandas as pd
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

def calculate_dinucleotide_metrics(file_path):
    """
    Calculate overall counts, rates, and O/E ratios for A, T, G, C and all 16 dinucleotides 
    for each gene region ("upstream1", "upstream2", "exon", "intron", "downstream1", "downstream2") in a single file.
    """
    try:
        gene_regions = ["upstream1", "upstream2", "exon", "intron", "downstream1", "downstream2"]
        bases = ['A', 'T', 'G', 'C']
        dinucleotides = [a + b for a in bases for b in bases]

        cumulative_metrics = {region: {"A": 0, "T": 0, "G": 0, "C": 0, "valid_dinucleotide_counts": 0, **{di: 0 for di in dinucleotides}} for region in gene_regions}
        
        df = pd.read_csv(file_path, sep="\t")
        df["Region"] = df["Busco ID"].apply(lambda x: x.split("_")[-1])

        for _, row in df.iterrows():
            region = row["Region"]
            if region in gene_regions:
                cumulative_metrics[region]["A"] += row["A Count"]
                cumulative_metrics[region]["T"] += row["T Count"]
                cumulative_metrics[region]["G"] += row["G Count"]
                cumulative_metrics[region]["C"] += row["C Count"]
                cumulative_metrics[region]["valid_dinucleotide_counts"] += row["Valid Dinucleotide Counts"]
                for dinucleotide in dinucleotides:
                    cumulative_metrics[region][dinucleotide] += row.get(f"{dinucleotide} Count", 0)

        results = []
        for region, metrics in cumulative_metrics.items():
            total_length = metrics["valid_dinucleotide_counts"]
            if total_length > 0:
                total_base_count = metrics["A"] + metrics["T"] + metrics["G"] + metrics["C"]
                if total_base_count > 0:
                    base_freq = {
                        nuc: metrics[nuc] / total_base_count for nuc in ["A", "T", "G", "C"]
                    }
                    rates = base_freq.copy()
                    expected_counts = {
                        di: base_freq[di[0]] * base_freq[di[1]] * total_length
                        for di in dinucleotides
                    }
                else:
                    base_freq = {nuc: 0 for nuc in ["A", "T", "G", "C"]}
                    rates = base_freq.copy()
                    expected_counts = {di: 0 for di in dinucleotides}

                # Compute O/E ratios
                oe_ratios = {
                    di: (metrics[di] / expected_counts[di]) if expected_counts[di] > 0 else 0
                    for di in dinucleotides
                }

                dinucleotide_rates = {di: metrics[di] / total_length for di in dinucleotides}

                results.append({
                    "Species": os.path.basename(file_path).replace("_gene_summary.txt", ""),
                    "Region": region,
                    "Total Valid Dinucleotide Counts": total_length,
                    "A Count": metrics["A"],
                    "T Count": metrics["T"],
                    "G Count": metrics["G"],
                    "C Count": metrics["C"],
                    "A Rate": rates["A"],
                    "T Rate": rates["T"],
                    "G Rate": rates["G"],
                    "C Rate": rates["C"],
                    **{f"{di} Count": metrics[di] for di in dinucleotides},
                    **{f"{di} Rate": dinucleotide_rates[di] for di in dinucleotides},
                    **{f"{di} O/E": oe_ratios[di] for di in dinucleotides},
                })
        return results
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return []

def summarize_all_species(input_dir, output_file, num_processes=None):
    """
    Summarize metrics for all species into a single file using multiprocessing.
    """
    all_results = []
    files = [os.path.join(input_dir, file) for file in os.listdir(input_dir) if file.endswith("_gene_summary.txt")]
    
    if num_processes is None:
        num_processes = max(1, os.cpu_count() - 1)

    print(f"Processing {len(files)} files using {num_processes} parallel processes...")

    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        future_to_file = {executor.submit(calculate_dinucleotide_metrics, file): file for file in files}

        for future in as_completed(future_to_file):
            file = future_to_file[future]
            try:
                results = future.result()
                all_results.extend(results)
                print(f"Completed: {os.path.basename(file)}")
            except Exception as e:
                print(f"Error processing {file}: {e}")

    all_results_df = pd.DataFrame(all_results)
    all_results_df.to_csv(output_file, sep=",", index=False)
    print(f"Combined results saved to {output_file}")

if __name__ == "__main__":
    input_dir = "/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/dinucleotide_count/busco_gene_summaries_Github_check"
    output_file = "/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/dinucleotide_count/summary_dinucleotide_rates_region_separated_Github_check.csv"

    num_processes = max(1, os.cpu_count() - 1)
    summarize_all_species(input_dir, output_file, num_processes)

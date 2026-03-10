import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed

REGION_LABEL = "gene_associated_-2kb_+1kb"
BASES = ["A", "T", "G", "C"]
DINUCS = [a + b for a in BASES for b in BASES]

def summarize_one_species(file_path: str):
    try:
        df = pd.read_csv(file_path, sep="\t")

        # Required columns in your actual files
        required = ["A_Count", "T_Count", "G_Count", "C_Count", "Valid_Dinucleotide_Counts"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns {missing}. Found: {list(df.columns)}")

        # Sum mono counts and valid dinuc denominator
        A = df["A_Count"].sum()
        T = df["T_Count"].sum()
        G = df["G_Count"].sum()
        C = df["C_Count"].sum()
        total_valid_dinuc = df["Valid_Dinucleotide_Counts"].sum()

        total_base = A + T + G + C
        if total_base > 0:
            mono_freq = {"A": A/total_base, "T": T/total_base, "G": G/total_base, "C": C/total_base}
        else:
            mono_freq = {b: 0 for b in BASES}

        # Sum dinucleotide counts
        dinuc_counts = {}
        for di in DINUCS:
            col = f"{di}_Count"
            dinuc_counts[di] = df[col].sum() if col in df.columns else 0

        # Dinucleotide rates (observed / total valid dinuc)
        dinuc_rates = {di: (dinuc_counts[di] / total_valid_dinuc) if total_valid_dinuc > 0 else 0 for di in DINUCS}

        # Expected counts from summed mono freqs
        expected = {di: (mono_freq[di[0]] * mono_freq[di[1]] * total_valid_dinuc) if total_valid_dinuc > 0 else 0 for di in DINUCS}
        dinuc_oe = {di: (dinuc_counts[di] / expected[di]) if expected[di] > 0 else 0 for di in DINUCS}

        # Species name from filename
        species = os.path.basename(file_path).replace("_genomic_gene_summary.txt", "")

        row = {
            "Species": species,
            "Region": REGION_LABEL,
            "Total_Valid_Dinucleotide_Counts": total_valid_dinuc,
            "A_Count": A, "T_Count": T, "G_Count": G, "C_Count": C,
            "A_Rate": mono_freq["A"], "T_Rate": mono_freq["T"], "G_Rate": mono_freq["G"], "C_Rate": mono_freq["C"],
        }

        for di in DINUCS:
            row[f"{di}_Count"] = dinuc_counts[di]
            row[f"{di}_Rate"] = dinuc_rates[di]
            row[f"{di}_OE"] = dinuc_oe[di]

        return row

    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return None

def summarize_all_species(input_dir: str, output_file: str, num_processes=None):
    files = [
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.endswith("_genomic_gene_summary.txt")
    ]

    if num_processes is None:
        num_processes = max(1, os.cpu_count() - 1)

    print(f"Processing {len(files)} files using {num_processes} processes...")

    rows = []
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        future_to_file = {executor.submit(summarize_one_species, fp): fp for fp in files}
        for future in as_completed(future_to_file):
            fp = future_to_file[future]
            res = future.result()
            if res is not None:
                rows.append(res)
            print(f"Completed: {os.path.basename(fp)}")

    out_df = pd.DataFrame(rows)
    out_df.to_csv(output_file, index=False)
    print(f"Saved: {output_file} (rows={len(out_df)})")

if __name__ == "__main__":
    input_dir = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/dinucleotide_count_merged/augustus_gene_summaries_Github_check"
    output_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/dinucleotide_count_merged/summary_dinucleotide_rates_gene_associated_-2kb_+1kb.csv"

    summarize_all_species(input_dir, output_file, num_processes=max(1, os.cpu_count() - 1))

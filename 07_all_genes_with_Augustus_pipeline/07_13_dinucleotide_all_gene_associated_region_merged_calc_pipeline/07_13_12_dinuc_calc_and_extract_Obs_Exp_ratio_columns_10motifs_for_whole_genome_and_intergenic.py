import pandas as pd
import numpy as np

# Load the input CSV file
input_file = "filtered_species_with_zero_lifecycle_6_categories_Github_check_with_intergenic_counts_rates_OE_with_tip_name.csv"
df = pd.read_csv(input_file)

# ---- Helper: safely get a numeric Series (or None if missing) ----
def get_numeric(series_name: str):
    if series_name not in df.columns:
        return None
    return pd.to_numeric(df[series_name], errors="coerce")

# ---- Define which rows to collapse (NOW: all regions) ----
# i.e., collapse reverse complements for genome, gene_associated_-2kb_+1kb, and intergenic
target_regions = {"genome", "gene_associated_-2kb_+1kb", "intergenic"}
is_target = df["Region"].astype(str).isin(target_regions)

# ---- Define reverse-complement groups for dinucleotides ----
# Pairs: collapse with formula using observed frequencies (Rates)
#   OE_collapsed = (f_XY + f_RC) / (f_X f_Y + f_R f_C)
# Palindromes: AT, TA, CG, GC are their own reverse complements (keep as-is)
rc_pairs = [
    # (dinuc1, dinuc2, new_name)
    ("TG", "CA", "CA_and_TG_OE"),
    ("AC", "GT", "AC_and_GT_OE"),
    ("AG", "CT", "AG_and_CT_OE"),
    ("GA", "TC", "GA_and_TC_OE"),
    ("AA", "TT", "AA_and_TT_OE"),
    ("CC", "GG", "CC_and_GG_OE"),
]
palindromes = ["AT", "TA", "CG", "GC"]

# ---- Compute collapsed O/E metrics for target rows ----
# Requires: mononucleotide frequencies A_Rate,C_Rate,G_Rate,T_Rate and dinucleotide frequencies XY_Rate
A = get_numeric("A_Rate")
C = get_numeric("C_Rate")
G = get_numeric("G_Rate")
T = get_numeric("T_Rate")

have_mono = all(x is not None for x in [A, C, G, T])

if not have_mono:
    print("WARNING: Missing one or more mononucleotide frequency columns: A_Rate, C_Rate, G_Rate, T_Rate.")
    print("         Collapsed O/E columns will not be computed.")
else:
    mono = {"A": A, "C": C, "G": G, "T": T}

    for d1, d2, newcol in rc_pairs:
        f1 = get_numeric(f"{d1}_Rate")
        f2 = get_numeric(f"{d2}_Rate")
        if f1 is None or f2 is None:
            print(f"WARNING: Missing {d1}_Rate or {d2}_Rate. Skipping {newcol}.")
            continue

        # Bases for denominator terms:
        # d1 = b1b2, d2 = r1r2
        b1, b2 = d1[0], d1[1]
        r1, r2 = d2[0], d2[1]

        denom = (mono[b1] * mono[b2]) + (mono[r1] * mono[r2])
        numer = f1 + f2

        collapsed = numer / denom

        # Only fill target rows; other regions stay NaN
        if newcol not in df.columns:
            df[newcol] = np.nan
        df.loc[is_target, newcol] = collapsed.loc[is_target]

        # Blank out the original paired OE columns for target rows
        # so those regions effectively use 10 collapsed metrics
        oe1 = f"{d1}_OE"
        oe2 = f"{d2}_OE"
        if oe1 in df.columns:
            df.loc[is_target, oe1] = np.nan
        if oe2 in df.columns:
            df.loc[is_target, oe2] = np.nan

    # Palindromes: do nothing (keep existing AT_OE, TA_OE, CG_OE, GC_OE)

# ---- Now filter to output columns: metadata + *_OE (including new collapsed ones) ----
meta_cols = ["Organism Name", "Region", "Average_lifespan_days", "Tip_Name", "Genome_size_bp"]

# Keep only columns that exist (avoid KeyError if any meta col missing)
meta_cols_existing = [c for c in meta_cols if c in df.columns]

columns_to_keep = [
    col for col in df.columns
    if ("_OE" in col) or (col in meta_cols_existing)
]
df_filtered = df[columns_to_keep].copy()

# Remove mono-nucleotide OE columns if present
columns_to_remove = ["A_OE", "T_OE", "G_OE", "C_OE"]
df_filtered = df_filtered.drop(
    columns=[c for c in columns_to_remove if c in df_filtered.columns],
    errors="ignore"
)

# Reorder columns: metadata first
ordered_columns = meta_cols_existing + [c for c in df_filtered.columns if c not in meta_cols_existing]
df_filtered = df_filtered[ordered_columns]

# Save
output_file = "di_obs_exp_ratio_rev_comp_collapsed_only_with_lifespan_tip_names_and_genome_size_Github_check.csv"
df_filtered.to_csv(output_file, index=False)

print(f"Filtered data saved to {output_file}")
print(f"Collapsed reverse-complement dinucleotide O/E for regions: {sorted(target_regions)}")

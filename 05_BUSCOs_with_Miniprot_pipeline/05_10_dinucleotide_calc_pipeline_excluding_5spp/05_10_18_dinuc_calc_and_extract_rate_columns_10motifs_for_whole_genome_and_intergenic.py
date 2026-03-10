import pandas as pd
import numpy as np

# Load the input CSV file
input_file = "filtered_species_with_tip_name_dinuc_Github_check_plus_intergenic.csv"
df = pd.read_csv(input_file)

# Retain only the columns you want
meta_cols = ["Organism Name", "Region", "Average_lifespan_days", "Tip_Name", "Genome_size_bp"]
columns_to_keep = [c for c in df.columns if (c in meta_cols) or ("_Rate" in c)]
df_filtered = df[columns_to_keep].copy()

# Remove mono-nucleotide rate columns if present
columns_to_remove = ["A_Rate", "T_Rate", "G_Rate", "C_Rate"]
df_filtered = df_filtered.drop(
    columns=[c for c in columns_to_remove if c in df_filtered.columns],
    errors="ignore"
)

# --- Reverse-complement collapsing for Region in {"genome", "intergenic"} ---
# 16 dinucs -> 10 RC groups:
#   6 paired groups become new columns:
#     {AA,TT}, {AC,GT}, {AG,CT}, {CA,TG}, {CC,GG}, {GA,TC}
#   4 palindromic motifs are kept AS-IS (no renaming):
#     AT_Rate, TA_Rate, CG_Rate, GC_Rate
rc_pairs = [
    ("AA_Rate", "TT_Rate", "AA_and_TT_Rate"),
    ("AC_Rate", "GT_Rate", "AC_and_GT_Rate"),
    ("AG_Rate", "CT_Rate", "AG_and_CT_Rate"),
    ("CA_Rate", "TG_Rate", "CA_and_TG_Rate"),
    ("CC_Rate", "GG_Rate", "CC_and_GG_Rate"),
    ("GA_Rate", "TC_Rate", "GA_and_TC_Rate"),
]

palindromes_keep = ["AT_Rate", "TA_Rate", "CG_Rate", "GC_Rate"]

is_target = df_filtered["Region"].astype(str).isin(["genome", "intergenic"])

# Create collapsed columns for target rows; blank out originals for those pairs in target rows
for col1, col2, newcol in rc_pairs:
    if col1 not in df_filtered.columns or col2 not in df_filtered.columns:
        continue

    if newcol not in df_filtered.columns:
        df_filtered[newcol] = np.nan

    df_filtered.loc[is_target, newcol] = (
        pd.to_numeric(df_filtered.loc[is_target, col1], errors="coerce") +
        pd.to_numeric(df_filtered.loc[is_target, col2], errors="coerce")
    )

    # For target rows, blank out the original paired columns so those regions effectively use 10 motifs
    df_filtered.loc[is_target, [col1, col2]] = np.nan

# IMPORTANT:
# - We do NOT create AT_Rate_rc / TA_Rate_rc / CG_Rate_rc / GC_Rate_rc.
# - Palindromic columns remain as AT_Rate, TA_Rate, CG_Rate, GC_Rate.

# Reorder columns: metadata first, then the rest
ordered_columns = meta_cols + [c for c in df_filtered.columns if c not in meta_cols]
df_filtered = df_filtered[ordered_columns]

# Save the filtered data to a new CSV file
output_file = "di_rates_rev_comp_collapsed_genome_and_intergenic_with_lifespan_tip_names_and_genome_size_Github_check.csv"
df_filtered.to_csv(output_file, index=False)

print(f"Filtered data saved to {output_file}")

import pandas as pd
import numpy as np

# ---- INPUT / OUTPUT ----
input_file  = "filtered_species_with_zero_lifecycle_6_categories_Github_check.csv"
output_file = "filtered_species_with_zero_lifecycle_6_categories_Github_check_with_intergenic_counts_rates_OE.csv"

# ---- LOAD ----
df = pd.read_csv(input_file)

# Species key (prefer Assembly Identifier if present)
species_key = "Assembly Identifier" if "Assembly Identifier" in df.columns else "Organism Name"
if species_key not in df.columns:
    raise ValueError('Neither "Assembly Identifier" nor "Organism Name" exists in the CSV.')

# ---- REGION LABELS ----
REGION_GENOME = "genome"
REGION_GENE   = "gene_associated_-2kb_+1kb"
REGION_INTER  = "intergenic"

# ---- BASES / COLUMNS ----
bases = ["A", "C", "G", "T"]

# Mononucleotide count/rate columns
mono_count_cols = [f"{b}_Count" for b in bases]
missing_mono_counts = [c for c in mono_count_cols if c not in df.columns]
if missing_mono_counts:
    raise ValueError(f"Missing mononucleotide count columns: {missing_mono_counts}")

mono_rate_cols = [f"{b}_Rate" for b in bases]

# Dinucleotide count/rate/OE columns
dinuc_pairs = [b1 + b2 for b1 in bases for b2 in bases]
dinuc_count_cols = [f"{p}_Count" for p in dinuc_pairs]
missing_dinuc_counts = [c for c in dinuc_count_cols if c not in df.columns]
if missing_dinuc_counts:
    raise ValueError(f"Missing dinucleotide count columns (examples): {missing_dinuc_counts[:8]}")

dinuc_rate_cols = [f"{p}_Rate" for p in dinuc_pairs]
dinuc_oe_cols   = [f"{p}_OE"   for p in dinuc_pairs]

# ---- COERCE COUNTS TO NUMERIC ----
df[mono_count_cols + dinuc_count_cols] = df[mono_count_cols + dinuc_count_cols].apply(
    pd.to_numeric, errors="coerce"
)

# ---- MAKE INTERGENIC ROWS (COUNTS) ----
count_cols = mono_count_cols + dinuc_count_cols
new_rows = []

for sp, g in df.groupby(species_key, dropna=False):
    genome_row = g[g["Region"] == REGION_GENOME]
    gene_row   = g[g["Region"] == REGION_GENE]

    if len(genome_row) != 1 or len(gene_row) != 1:
        continue

    genome_row = genome_row.iloc[0]
    gene_row   = gene_row.iloc[0]

    intergenic_counts = (genome_row[count_cols] - gene_row[count_cols]).clip(lower=0)

    new_row = genome_row.copy()  # copy metadata from genome row
    new_row["Region"] = REGION_INTER
    new_row[count_cols] = intergenic_counts

    new_rows.append(new_row)

df_out = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)

# ---- RETAIN EXISTING Sum_of_valid_dinucleotide_counts ----
if "Sum_of_valid_dinucleotide_counts" not in df_out.columns:
    df_out["Sum_of_valid_dinucleotide_counts"] = np.nan

# ---- ENSURE RATE + OE COLUMNS EXIST (DON'T OVERWRITE NON-INTERGENIC) ----
for c in mono_rate_cols + dinuc_rate_cols + dinuc_oe_cols:
    if c not in df_out.columns:
        df_out[c] = np.nan

# ---- CALCULATE RATES + SUM ONLY FOR INTERGENIC ROWS ----
m_inter = df_out["Region"] == REGION_INTER

# Mononucleotide rates
mono_sum = df_out.loc[m_inter, mono_count_cols].sum(axis=1).replace({0: np.nan})
for b in bases:
    df_out.loc[m_inter, f"{b}_Rate"] = df_out.loc[m_inter, f"{b}_Count"] / mono_sum

# Dinucleotide rates + sum of valid dinucleotide counts
dinuc_sum = df_out.loc[m_inter, dinuc_count_cols].sum(axis=1)
dinuc_sum_safe = dinuc_sum.replace({0: np.nan})

df_out.loc[m_inter, "Sum_of_valid_dinucleotide_counts"] = dinuc_sum

for p in dinuc_pairs:
    df_out.loc[m_inter, f"{p}_Rate"] = df_out.loc[m_inter, f"{p}_Count"] / dinuc_sum_safe

# ---- CALCULATE DINUCLEOTIDE O/E ONLY FOR INTERGENIC ROWS ----
# OE(XY) = XY_Rate / (X_Rate * Y_Rate)
for p in dinuc_pairs:
    x, y = p[0], p[1]
    denom = df_out.loc[m_inter, f"{x}_Rate"] * df_out.loc[m_inter, f"{y}_Rate"]
    denom = denom.replace({0: np.nan})
    df_out.loc[m_inter, f"{p}_OE"] = df_out.loc[m_inter, f"{p}_Rate"] / denom

# ---- SAVE ----
df_out.to_csv(output_file, index=False)

print("Done.")
print(f"  Added {len(new_rows)} intergenic rows (Region='{REGION_INTER}').")
print("  Filled mono/dinuc rates + dinuc O/E for intergenic rows only.")
print("  Retained existing values for non-intergenic rows (including Sum_of_valid_dinucleotide_counts).")
print(f"Wrote: {output_file}")

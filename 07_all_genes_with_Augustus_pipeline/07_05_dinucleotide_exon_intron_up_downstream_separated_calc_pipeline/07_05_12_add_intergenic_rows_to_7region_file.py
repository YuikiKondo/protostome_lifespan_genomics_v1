import pandas as pd

# ---- INPUTS ----
augustus_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/gene_associated_-2kb_+1kb/dinucleotide_count_merged/dinucleotide_PGLS/filtered_species_with_zero_lifecycle_6_categories_Github_check_with_intergenic_counts_rates_OE_with_gene_associated_total_bp.csv" #output file from 07_13_14_add_total_gene_region_bp_count.py

busco_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/up_downstreams_range_refined_exon_intron_separated/dinucleotide_count/dinucleotide_PGLS/filtered_species_with_tip_name_dinuc_Github_check.csv"

# ---- OUTPUT ----
output_file = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/up_downstreams_range_refined_exon_intron_separated/dinucleotide_count/dinucleotide_PGLS/filtered_species_with_tip_name_dinuc_Github_check_plus_intergenic.csv"

# ---- LOAD ----
df_aug = pd.read_csv(augustus_file)
df_bus = pd.read_csv(busco_file)

# ---- REQUIRED COLUMN NAMES ----
new_col = "Total_bp_gene_associated_region_including_N"
genome_col = "Genome_size_bp"
asm_col = "Assembly Identifier"
region_col = "Region"

# ---- sanity checks ----
required_aug = {asm_col, region_col, new_col}
missing_aug = required_aug - set(df_aug.columns)
if missing_aug:
    raise ValueError(f"Augustus file missing required columns: {missing_aug}")

required_bus = {asm_col, region_col, genome_col}
missing_bus = required_bus - set(df_bus.columns)
if missing_bus:
    raise ValueError(f"BUSCO file missing required columns: {missing_bus}")

# ---- Extract intergenic rows from Augustus ----
intergenic = df_aug[df_aug[region_col].astype(str) == "intergenic"].copy()
if intergenic.empty:
    raise ValueError(f"No intergenic rows found in: {augustus_file}")

# ---- Ensure BUSCO has new columns so output schema includes them ----
if new_col not in df_bus.columns:
    df_bus[new_col] = pd.NA

intergenic_bp_col = "Intergenic_bp"
if intergenic_bp_col not in df_bus.columns:
    df_bus[intergenic_bp_col] = pd.NA

# ---- ALIGN intergenic columns to BUSCO schema (robust) ----
bus_cols = list(df_bus.columns)

for c in bus_cols:
    if c not in intergenic.columns:
        intergenic[c] = pd.NA

extra_in_intergenic = [c for c in intergenic.columns if c not in bus_cols]
if extra_in_intergenic:
    intergenic = intergenic.drop(columns=extra_in_intergenic)

intergenic = intergenic[bus_cols]

# ---- OPTIONAL: avoid duplicates if BUSCO already has intergenic ----
mask_existing_intergenic = df_bus[region_col].astype(str) == "intergenic"
if mask_existing_intergenic.any():
    key_cols = [asm_col, region_col]
    incoming_keys = set(tuple(x) for x in intergenic.loc[:, key_cols].astype(str).values)

    to_drop = []
    for idx, row in df_bus.loc[mask_existing_intergenic, key_cols].astype(str).iterrows():
        if tuple(row.values) in incoming_keys:
            to_drop.append(idx)

    if to_drop:
        df_bus = df_bus.drop(index=to_drop)

# ---- CONCATENATE ----
df_out = pd.concat([df_bus, intergenic], ignore_index=True)

# ---- Broadcast Total_bp_gene_associated_region_including_N to ALL regions per species ----
# Build mapping from Augustus file (unique per Assembly Identifier)
map_df = (
    df_aug[[asm_col, new_col]]
    .dropna(subset=[asm_col, new_col])
    .drop_duplicates(subset=[asm_col])
)

bp_map = dict(zip(map_df[asm_col].astype(str), map_df[new_col]))

# Fill / overwrite df_out[new_col] using mapping (for ALL rows)
df_out[new_col] = df_out[asm_col].astype(str).map(bp_map)

# ---- Compute Intergenic_bp = Genome_size_bp - Total_bp_gene_associated... ----
df_out[genome_col] = pd.to_numeric(df_out[genome_col], errors="coerce")
df_out[new_col] = pd.to_numeric(df_out[new_col], errors="coerce")

df_out[intergenic_bp_col] = df_out[genome_col] - df_out[new_col]

# ---- Reorder columns: place new_col right next to Genome_size_bp, and Intergenic_bp right after ----
cols = list(df_out.columns)

# Remove if already present (we'll insert at desired position)
for c in [new_col, intergenic_bp_col]:
    if c in cols:
        cols.remove(c)

if genome_col in cols:
    gidx = cols.index(genome_col) + 1
    cols[gidx:gidx] = [new_col, intergenic_bp_col]
else:
    # fallback: append at end (shouldn't happen due to earlier check)
    cols.extend([new_col, intergenic_bp_col])

df_out = df_out[cols]

# ---- (Optional) sort for readability ----
sort_cols = [c for c in ["Tip_Name", "Organism Name", asm_col, region_col] if c in df_out.columns]
if sort_cols:
    df_out = df_out.sort_values(sort_cols, kind="mergesort").reset_index(drop=True)

# ---- Report any missing mapping (should be rare, but good to know) ----
n_missing_bp = df_out[new_col].isna().sum()
if n_missing_bp > 0:
    print(f"WARNING: {n_missing_bp} rows have missing {new_col} (no match by {asm_col}).")

# ---- SAVE ----
df_out.to_csv(output_file, index=False)

print("Done.")
print(f"  BUSCO rows (after any de-dup): {len(df_bus)}")
print(f"  Intergenic added:              {len(intergenic)}")
print(f"  Output rows total:             {len(df_out)}")
print(f"Wrote: {output_file}")
print(f"Added columns next to {genome_col}: {new_col}, {intergenic_bp_col}")

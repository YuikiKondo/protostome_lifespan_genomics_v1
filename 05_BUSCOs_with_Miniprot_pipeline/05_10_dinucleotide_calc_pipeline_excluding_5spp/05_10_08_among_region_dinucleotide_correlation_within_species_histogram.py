import pandas as pd
import matplotlib.pyplot as plt

# Load the results file
file_path = '/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_range_refined_exon_intron_separated_including_dup/gene_region_separated_sequences_all_hit_buscos/dinucleotide_count_excluding_5spp/among_region_dinucleotide_correlation_results_Github_check.csv'
correlation_data = pd.read_csv(file_path)

# Define regions and colors
regions = ["intergenic", "upstream1", "upstream2", "exon", "intron", "downstream1", "downstream2"]
region_colors = {
    "intergenic": "grey",
    "upstream1": "#4682B4",  # Steel Blue
    "upstream2": "#5F9EA0",  # Cadet Blue
    "exon": "#228B22",       # Forest Green
    "intron": "#FFA500",     # Orange
    "downstream1": "#DC143C", # Crimson
    "downstream2": "#B22222"  # Firebrick
}

# Set up the figure and subplots
fig, axes = plt.subplots(len(regions), 1, figsize=(10, 14), sharex=True)

# Plot the histograms
for i, region in enumerate(regions):
    region_data = correlation_data[correlation_data["Comparison"] == f"genome vs {region}"]
    total_species = len(region_data)
    bonferroni_significant = (region_data["Adjusted p-value (Bonferroni)"] < 0.05).sum()
    
    axes[i].hist(
        region_data["Correlation (r)"], 
        bins=30, 
        alpha=0.7, 
        color=region_colors[region]
    )
    
    # Add title with significant species count
    axes[i].set_title(
        f"Genome vs {region} (Bonferroni p<0.05: {bonferroni_significant} species out of {total_species} species)", 
        fontsize=12
    )
    axes[i].set_ylabel("Number of species")
    axes[i].grid(axis='y', alpha=0.75)

# Add common x-axis label
axes[-1].set_xlabel("Correlation r")

# Save the plot to the current working directory
output_path = "among_region_dinucleotide_correlation_histogram_Github_check.png"
plt.tight_layout()
plt.savefig(output_path)

print(f"Histogram saved as: {output_path}")

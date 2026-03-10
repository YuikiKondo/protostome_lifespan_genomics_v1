import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

# Load the dataset
file_path = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/up_downstreams_range_refined_exon_intron_separated/dinucleotide_count/filtered_dinucleotide_lifespan_class_Github_check.csv"
data = pd.read_csv(file_path)

# Filter data for organisms with lifespan data
filtered_data = data[data['Average_lifespan_days'].notnull()]

# Define the regions and dinucleotide rate columns
regions = ["genome", "upstream1", "upstream2", "exon", "intron", "downstream1", "downstream2"]
dinucleotide_rates = [
    "AA_Rate", "AT_Rate", "AG_Rate", "AC_Rate",
    "TA_Rate", "TT_Rate", "TG_Rate", "TC_Rate",
    "GA_Rate", "GT_Rate", "GG_Rate", "GC_Rate",
    "CA_Rate", "CT_Rate", "CG_Rate", "CC_Rate"
]

# Initialize a list to store results
correlation_results = []

# Calculate total number of tests for Bonferroni adjustment
num_species = filtered_data["Organism Name"].nunique()
total_tests = len(regions) * num_species  # 4 comparisons per species

# Iterate through each organism
for organism in filtered_data["Organism Name"].unique():
    organism_data = filtered_data[filtered_data["Organism Name"] == organism]
    genome_data = organism_data[organism_data["Region"] == "genome"]
    
    # Skip if no genome data available
    if genome_data.empty:
        continue

    genome_rates = genome_data[dinucleotide_rates].iloc[0]
    p_values = []
    comparisons = []
    r_values = []

    for region in regions:
        if region == "genome":
            continue  # Skip self-correlation
        
        region_data = organism_data[organism_data["Region"] == region]
        if region_data.empty:
            continue  # Skip if region data is unavailable

        region_rates = region_data[dinucleotide_rates].iloc[0]
        
        # Calculate correlation for each dinucleotide rate
        r, p = stats.pearsonr(genome_rates, region_rates)
        p_values.append(p)
        r_values.append(r)
        comparisons.append(f"genome vs {region}")
    
    # Adjust p-values for multiple comparisons
    if p_values:
        # Adjust using Bonferroni by the total number of tests
        bonferroni_adjusted_p = [min(p * total_tests, 1) for p in p_values]
        
        # Adjust using FDR (Benjamini-Hochberg) for all tests globally
        fdr_adjusted_p = multipletests(p_values, method="fdr_bh")[1]
        
        # Store the results
        for i, comparison in enumerate(comparisons):
            correlation_results.append({
                "Organism Name": organism,
                "Comparison": comparison,
                "Correlation (r)": r_values[i],
                "p-value": p_values[i],
                "Adjusted p-value (FDR)": fdr_adjusted_p[i],
                "Adjusted p-value (Bonferroni)": bonferroni_adjusted_p[i]
            })

# Convert results to a DataFrame
correlation_df = pd.DataFrame(correlation_results)

# Save to CSV
output_path = "/home/y-kondo/protostome_lifespan_model_augustus/Augustus_seq_directory/up_downstreams_range_refined_exon_intron_separated/dinucleotide_count/among_region_dinucleotide_correlation_results_Github_check.csv"
correlation_df.to_csv(output_path, index=False)

print(f"Results saved to: {output_path}")
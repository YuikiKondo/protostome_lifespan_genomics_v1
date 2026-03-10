# Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)

# Define file paths
input_file <- "/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/Downstream_analysis_Github_check/all_species_sliding_window_downstream_dinucleotide_counts_stepsize_10_with_species.csv"
output_dir <- "/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/Downstream_analysis_Github_check/Dinucleotide_rate_dinuc_plots_stepsize_10"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load the data
df <- read_csv(input_file)

# Define custom species order with new names
custom_species_names <- c(
  "Armandia brevis (Annelida short-lived)", "Hirudinaria manillensis (Annelida long-lived)",
  "Fuscozetes fuscipes (Chelicerata short-lived)", "Limulus polyphemus (Chelicerata long-lived)",
  "Anthophora plumipes (Mandibulata short-lived)", "Birgus latro (Mandibulata long-lived)",
  "Bullina lineata (Mollusca short-lived)", "Margaritifera margaritifera (Mollusca long-lived)",
  "Auanema rhodensis (Nematoda short-lived)", "Wuchereria bancrofti (Nematoda long-lived)",
  "Gyrodactylus salaris (Platyhelminthes short-lived)", "Taenia asiatica (Platyhelminthes long-lived)"
)

# Define a mapping from original species names to new names
species_mapping <- setNames(custom_species_names, unique(df$Species))

# Replace Species column with the new names
df <- df %>% mutate(Species = species_mapping[Species])

# Ensure factor levels for correct ordering
df$Species <- factor(df$Species, levels = custom_species_names)

# Get unique assemblies in the custom order
assemblies <- df %>%
  select(Assembly, Species) %>%
  distinct() %>%
  arrange(factor(Species, levels = custom_species_names)) %>%
  pull(Assembly)

# Define dinucleotide rate indices
dinucleotide_rates <- c("AA_rate", "AT_rate", "AG_rate", "AC_rate",
"TA_rate", "TT_rate", "TG_rate", "TC_rate",
"GA_rate", "GT_rate", "GG_rate", "GC_rate",
"CA_rate", "CT_rate", "CG_rate", "CC_rate")

# Generate one PNG file per DNA mechanics index
for (index in dinucleotide_rates) {
  
  # Calculate global y-axis limits using full min and max values
  ylim_min <- min(df[[index]], na.rm = TRUE)
  ylim_max <- max(df[[index]], na.rm = TRUE)

  # Store all species plots for this index
  plot_list <- list()
  
  for (assembly_name in assemblies) {
    subset_df <- df %>% filter(Assembly == assembly_name)
    
    # Extract species name corresponding to the assembly
    species_name <- unique(subset_df$Species)
    
    # Create the plot for the current index
    p <- ggplot(subset_df, aes(x = Position, y = .data[[index]])) +
      geom_line(color = "blue") +
      labs(
        x = "Distance from start codon (base)",
        y = paste(index, "index"),
        title = species_name
      ) +
      theme_classic() +
      coord_cartesian(ylim = c(ylim_min, ylim_max)) +  # Use full global limits
      geom_vline(xintercept = -100, lty = 2, colour = "grey10", linewidth = 0.5) +
      geom_vline(xintercept = 0, lty = 2, colour = "grey10", linewidth = 0.5) +
      geom_vline(xintercept = 200, lty = 2, colour = "grey10", linewidth = 0.5)
    
    plot_list[[species_name]] <- p
  }
  
  # Arrange all species in a grid (2 species per row)
  combined_plot <- ggarrange(plotlist = plot_list, ncol = 2, nrow = ceiling(length(plot_list) / 2))
  
  # Save each index plot as a separate PNG file
  output_file <- paste0(output_dir, "/", index, "_combined_plot.png")
  ggsave(output_file, combined_plot, width = 14, height = 6 * ceiling(length(plot_list) / 2), dpi = 300)
  
  print(paste("Saved:", output_file))
}

print("All dinucleotide rate index plots have been saved.")

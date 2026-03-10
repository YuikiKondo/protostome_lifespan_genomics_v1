## SET UP ----

## Load libraries
library(Biostrings)
library(tidyverse)
library(doParallel)
library(foreach)

## Define species list
species_list <- c(
  "GCA_035585475.1_ASM3558547v1",
  "GCA_034509925.1_ASM3450992v1",
  "GCA_034702105.1_ASM3470210v1",
  "GCA_000517525.1_Limulus_polyphemus-2.1.2",
  "GCA_951804975.1_iyAntPlum1.1",
  "GCA_018397915.1_ASM1839791v1",
  "GCA_039654405.1_ASM3965440v1",
  "GCA_029931535.1_MarmarV2",
  "GCA_964057225.1_nxAuaRhod1.1",
  "GCA_005281725.1_ASM528172v1",
  "GCA_000715275.1_Gsalaris_v1",
  "GCA_001693035.2_Taenia_asiatica_TASYD01_v1"
)

## Set window and step size for sliding window analysis
window_size <- 100
step_size <- 10

## Define all 16 dinucleotides
dinucleotides <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", 
                   "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")

## Parallel setup
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

## SLIDING WINDOW ANALYSIS ----

all_results <- foreach(assembly_name = species_list, .combine = bind_rows, .packages = c("Biostrings", "tidyverse")) %dopar% {
  species_files <- list(
    "complete" = paste0("/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/TSS_Github_check/representative_species_TSS_coordinates/", assembly_name, "_genomic.fna_complete_TSS_regions.fasta"),
    "duplicated" = paste0("/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/TSS_Github_check/representative_species_TSS_coordinates/", assembly_name, "_genomic.fna_duplicated_TSS_regions.fasta"),
    "fragmented" = paste0("/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/TSS_Github_check/representative_species_TSS_coordinates/", assembly_name, "_genomic.fna_fragmented_TSS_regions.fasta")
  )
  
  species_dinuc_results <- list()
  
  for (file_path in species_files) {
    if (!file.exists(file_path)) {
      message("File not found: ", file_path)
      next
    }
    
    temp_sequences <- readDNAStringSet(file_path)
    
    for (seq_idx in 1:length(temp_sequences)) {
      seq <- temp_sequences[[seq_idx]]
      
      if (length(seq) == 0) next
      
      seq_length <- length(seq)
      num_windows <- floor((seq_length - window_size) / step_size) + 1
      
      if (num_windows <= 0) next
      
      dinuc_counts_seq <- matrix(0, nrow = num_windows, ncol = length(dinucleotides))
      colnames(dinuc_counts_seq) <- paste0(dinucleotides, "_count")
      
      for (i in seq_len(num_windows)) {
        window_start <- (i - 1) * step_size + 1
        window_seq <- subseq(seq, start = window_start, end = window_start + window_size - 1)
        
        counts <- sapply(dinucleotides, function(dinuc) countPattern(dinuc, window_seq))
        dinuc_counts_seq[i, ] <- counts
      }
      
      df_counts <- data.frame(
        Assembly = assembly_name,
        Position = seq(-3000, by = step_size, length.out = num_windows),
        as.data.frame(dinuc_counts_seq)
      )
      
      df_counts$Total_Dinucleotide_Count <- rowSums(df_counts[, paste0(dinucleotides, "_count")], na.rm = TRUE)
      
      for (dinuc in dinucleotides) {
        df_counts[[paste0(dinuc, "_rate")]] <- df_counts[[paste0(dinuc, "_count")]] / df_counts$Total_Dinucleotide_Count
      }
      
      species_dinuc_results[[paste0(seq_idx, "_", basename(file_path))]] <- df_counts
    }
  }
  
  bind_rows(species_dinuc_results)
}

stopCluster(cl)

## DATA HANDLING ----

mean_dinuc_results <- all_results %>%
  group_by(Assembly, Position) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

## Save results
output_file <- "all_species_sliding_window_stepsize_10_dinucleotide_counts.csv"
write_csv(mean_dinuc_results, output_file)

print(paste("Results saved to", output_file))

## end script

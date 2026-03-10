## LOAD LIBRARIES ----
library(ape)
library(nlme)
library(tidyverse)

## LOAD DATA ----
input_file <- "filtered_species_with_tip_name_dinuc_Github_check_plus_intergenic.csv"
data <- read.csv(input_file, check.names = FALSE)

## BASIC PREP ----
data$Average_lifespan_days <- suppressWarnings(as.numeric(data$Average_lifespan_days))
data$ln_lifespan <- ifelse(!is.na(data$Average_lifespan_days) & data$Average_lifespan_days > 0,
                           log(data$Average_lifespan_days), NA_real_)

data$Genome_size_bp <- suppressWarnings(as.numeric(data$Genome_size_bp))
data$ln_Genome_Size <- ifelse(!is.na(data$Genome_size_bp) & data$Genome_size_bp > 0,
                              log(data$Genome_size_bp), NA_real_)

## Regions (8) ----
regions_to_keep <- c("genome", "intergenic", "upstream1", "upstream2",
                     "exon", "intron", "downstream1", "downstream2")
data$Region <- as.character(data$Region)
data <- data %>% filter(Region %in% regions_to_keep)

## GC content per row (region-specific, already normalized) ----
# Uses G_Rate + C_Rate
data$G_Rate <- suppressWarnings(as.numeric(data$G_Rate))
data$C_Rate <- suppressWarnings(as.numeric(data$C_Rate))
data$GC_content <- data$G_Rate + data$C_Rate

## Load tree ----
phylo_tree <- read.tree("/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_exon_intron_separated/gene_region_separated_sequences_no_internal_stop_3divisible/tetranucleotide_count/DNA_mechanics/PGLS/labelled_supertree_ottnames.tre")

## PREPROCESS DATA ----
matched_data <- data %>%
  filter(!is.na(Tip_Name), Tip_Name %in% phylo_tree$tip.label) %>%
  rename(species_name = Tip_Name)

phylo_tree <- compute.brlen(phylo_tree, method = "Grafen")
phylo_tree <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, matched_data$species_name))

regions <- unique(matched_data$Region[!is.na(matched_data$Region)])

## RUN PGLS (GC content as response) ----
pgls_results_all_regions <- list()

for (region in regions) {
  cat("Running PGLS (GC_content as response) for region:", region, "\n")

  region_data <- matched_data %>%
    filter(Region == region) %>%
    filter(species_name %in% phylo_tree$tip.label)

  if (nrow(region_data) < 3) {
    message("Skipping region ", region, " due to insufficient species (", nrow(region_data), " species).")
    next
  }

  ## Region-wise centering of ln genome size and ln lifespan
  region_data <- region_data %>%
    mutate(
      ln_Genome_Size = suppressWarnings(as.numeric(ln_Genome_Size)),
      ln_lifespan    = suppressWarnings(as.numeric(ln_lifespan)),
      GC_content     = suppressWarnings(as.numeric(GC_content)),
      ln_Genome_Size_c = as.numeric(scale(ln_Genome_Size, center = TRUE, scale = FALSE)),
      ln_lifespan_c    = as.numeric(scale(ln_lifespan, center = TRUE, scale = FALSE))
    )

  region_df <- region_data %>%
    transmute(
      species_name     = species_name,
      ln_lifespan_c    = as.numeric(ln_lifespan_c),
      ln_Genome_Size_c = as.numeric(ln_Genome_Size_c),
      GC_content       = as.numeric(GC_content)
    ) %>%
    drop_na(ln_lifespan_c, ln_Genome_Size_c, GC_content) %>%
    filter(is.finite(ln_lifespan_c), is.finite(ln_Genome_Size_c), is.finite(GC_content))

  if (nrow(region_df) < 3) {
    message("Skipping region ", region, " due to <3 complete cases (n=", nrow(region_df), ").")
    next
  }

  ## Model: GC_content as response
  formula_str <- "GC_content ~ ln_lifespan_c + ln_Genome_Size_c"

  tryCatch({
    pgls_model <- gls(
      as.formula(formula_str),
      data = region_df,
      correlation = corPagel(1, phy = phylo_tree, form = ~species_name),
      method = "ML"
    )

    model_summary <- summary(pgls_model)$tTable

    ## Extract lambda
    lambda_value <- tryCatch(
      as.numeric(coef(pgls_model$modelStruct$corStruct, unconstrained = FALSE)),
      error = function(e) NA_real_
    )

    ## Term-name safe extraction
    get_term <- function(term, col) {
      if (term %in% rownames(model_summary)) model_summary[term, col] else NA_real_
    }

    term_lifespan <- "ln_lifespan_c"
    term_genome   <- "ln_Genome_Size_c"

    pgls_results_all_regions[[paste(region, "GC_content", sep = "_")]] <- data.frame(
      Region = region,
      Trait = "GC_content",

      Intercept = get_term("(Intercept)", "Value"),
      Std.Error_Intercept = get_term("(Intercept)", "Std.Error"),
      t.value_Intercept = get_term("(Intercept)", "t-value"),
      p.value_Intercept = get_term("(Intercept)", "p-value"),

      Coefficient_ln_lifespan_c = get_term(term_lifespan, "Value"),
      Std.Error_ln_lifespan_c = get_term(term_lifespan, "Std.Error"),
      t.value_ln_lifespan_c = get_term(term_lifespan, "t-value"),
      p.value_ln_lifespan_c = get_term(term_lifespan, "p-value"),

      Coefficient_ln_Genome_Size_c = get_term(term_genome, "Value"),
      Std.Error_ln_Genome_Size_c = get_term(term_genome, "Std.Error"),
      t.value_ln_Genome_Size_c = get_term(term_genome, "t-value"),
      p.value_ln_Genome_Size_c = get_term(term_genome, "p-value"),

      Lambda = lambda_value,
      N_species = length(unique(region_df$species_name)),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message("Error in region ", region, ": ", e$message)
  })
}

## COMBINE RESULTS + FDR ----
if (length(pgls_results_all_regions) > 0) {
  pgls_results <- bind_rows(pgls_results_all_regions) %>%
    mutate(
      p.adjusted.fdr_Intercept        = p.adjust(p.value_Intercept, method = "BH"),
      p.adjusted.fdr_ln_lifespan_c    = p.adjust(p.value_ln_lifespan_c, method = "BH"),
      p.adjusted.fdr_ln_Genome_Size_c = p.adjust(p.value_ln_Genome_Size_c, method = "BH")
    )

  out_file <- "pgls_results_8_regions_GCcontent_LS_CENTERED_G_CENTERED_no_interaction_GC_as_response.csv"
  write.csv(pgls_results, out_file, row.names = FALSE)
  message("Saved: ", out_file)
} else {
  message("No PGLS results generated.")
}

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
regions_to_keep <- c("upstream1", "upstream2",
                     "exon", "intron", "downstream1", "downstream2")
data$Region <- as.character(data$Region)
data <- data %>% filter(Region %in% regions_to_keep)

## Traits to analyze (responses) ----
traits <- c("A_Rate", "T_Rate", "G_Rate", "C_Rate")

# Coerce trait columns to numeric (safe even if already numeric)
for (tr in traits) {
  if (tr %in% names(data)) data[[tr]] <- suppressWarnings(as.numeric(data[[tr]]))
}

## Load tree ----
phylo_tree <- read.tree("/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_exon_intron_separated/gene_region_separated_sequences_no_internal_stop_3divisible/tetranucleotide_count/DNA_mechanics/PGLS/labelled_supertree_ottnames.tre")

## PREPROCESS DATA ----
matched_data <- data %>%
  filter(!is.na(Tip_Name), Tip_Name %in% phylo_tree$tip.label) %>%
  rename(species_name = Tip_Name)

phylo_tree <- compute.brlen(phylo_tree, method = "Grafen")
phylo_tree <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, matched_data$species_name))

regions <- unique(matched_data$Region[!is.na(matched_data$Region)])

## RUN PGLS (mono-nucleotide rates as response) ----
pgls_results_all <- list()

for (region in regions) {
  cat("Region:", region, "\n")

  region_data <- matched_data %>%
    filter(Region == region) %>%
    filter(species_name %in% phylo_tree$tip.label)

  if (nrow(region_data) < 3) {
    message("Skipping region ", region, " due to insufficient species (", nrow(region_data), ").")
    next
  }

  ## Region-wise centering of ln genome size and ln lifespan
  region_data <- region_data %>%
    mutate(
      ln_Genome_Size = suppressWarnings(as.numeric(ln_Genome_Size)),
      ln_lifespan    = suppressWarnings(as.numeric(ln_lifespan)),
      ln_Genome_Size_c = as.numeric(scale(ln_Genome_Size, center = TRUE, scale = FALSE)),
      ln_lifespan_c    = as.numeric(scale(ln_lifespan, center = TRUE, scale = FALSE))
    )

  for (trait in traits) {
    if (!trait %in% names(region_data)) {
      message("Trait ", trait, " not found; skipping.")
      next
    }

    cat("  Running PGLS (", trait, " as response)\n", sep = "")

    region_df <- region_data %>%
      transmute(
        species_name     = species_name,
        ln_lifespan_c    = as.numeric(ln_lifespan_c),
        ln_Genome_Size_c = as.numeric(ln_Genome_Size_c),
        response         = suppressWarnings(as.numeric(.data[[trait]]))
      ) %>%
      drop_na(ln_lifespan_c, ln_Genome_Size_c, response) %>%
      filter(is.finite(ln_lifespan_c), is.finite(ln_Genome_Size_c), is.finite(response))

    if (nrow(region_df) < 3) {
      message("  Skipping ", region, " / ", trait, " due to <3 complete cases (n=", nrow(region_df), ").")
      next
    }

    ## Model: trait as response
    formula_str <- "response ~ ln_lifespan_c * ln_Genome_Size_c"

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
      term_int_1    <- "ln_lifespan_c:ln_Genome_Size_c"
      term_int_2    <- "ln_Genome_Size_c:ln_lifespan_c"
      term_int      <- if (term_int_1 %in% rownames(model_summary)) term_int_1 else term_int_2

      key <- paste(region, trait, sep = "_")
      pgls_results_all[[key]] <- data.frame(
        Region = region,
        Trait  = trait,

        Intercept = get_term("(Intercept)", "Value"),
        Std.Error_Intercept = get_term("(Intercept)", "Std.Error"),
        t.value_Intercept = get_term("(Intercept)", "t-value"),
        p.value_Intercept = get_term("(Intercept)", "p-value"),

        Coefficient_ln_lifespan_c = get_term(term_lifespan, "Value"),
        Std.Error_ln_lifespan_c   = get_term(term_lifespan, "Std.Error"),
        t.value_ln_lifespan_c     = get_term(term_lifespan, "t-value"),
        p.value_ln_lifespan_c     = get_term(term_lifespan, "p-value"),

        Coefficient_ln_Genome_Size_c = get_term(term_genome, "Value"),
        Std.Error_ln_Genome_Size_c   = get_term(term_genome, "Std.Error"),
        t.value_ln_Genome_Size_c     = get_term(term_genome, "t-value"),
        p.value_ln_Genome_Size_c     = get_term(term_genome, "p-value"),

        Coefficient_Interaction = get_term(term_int, "Value"),
        Std.Error_Interaction   = get_term(term_int, "Std.Error"),
        t.value_Interaction     = get_term(term_int, "t-value"),
        p.value_Interaction     = get_term(term_int, "p-value"),

        Lambda   = lambda_value,
        N_species = length(unique(region_df$species_name)),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      message("  Error in ", region, " / ", trait, ": ", e$message)
    })
  }
}

## COMBINE RESULTS + FDR ----
if (length(pgls_results_all) > 0) {
  pgls_results <- bind_rows(pgls_results_all) %>%
    mutate(
      p.adjusted.fdr_Intercept        = p.adjust(p.value_Intercept, method = "BH"),
      p.adjusted.fdr_ln_lifespan_c    = p.adjust(p.value_ln_lifespan_c, method = "BH"),
      p.adjusted.fdr_ln_Genome_Size_c = p.adjust(p.value_ln_Genome_Size_c, method = "BH"),
      p.adjusted.fdr_Interaction      = p.adjust(p.value_Interaction, method = "BH")
    )

  out_file <- "pgls_results_6_gene_associated_regions_mononuc_rates_LS_CENTERED_G_CENTERED_interaction_mono_as_response.csv"
  write.csv(pgls_results, out_file, row.names = FALSE)
  message("Saved: ", out_file)
} else {
  message("No PGLS results generated.")
}

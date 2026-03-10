## LOAD LIBRARIES ----
library(ape)
library(nlme)
library(tidyverse)

## LOAD DATA ----
input_file <- "di_obs_exp_ratio_rev_comp_collapsed_only_with_lifespan_tip_names_and_genome_size_Github_check.csv"
data <- read.csv(input_file, check.names = FALSE)

## BASIC PREP ----
data$ln_lifespan <- log(suppressWarnings(as.numeric(data$Average_lifespan_days)))

data$Genome_size_bp <- suppressWarnings(as.numeric(data$Genome_size_bp))
data$ln_Genome_Size <- ifelse(!is.na(data$Genome_size_bp) & data$Genome_size_bp > 0,
                              log(data$Genome_size_bp), NA_real_)

regions <- unique(data$Region[!is.na(data$Region)])

phylo_tree <- read.tree("/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_exon_intron_separated/gene_region_separated_sequences_no_internal_stop_3divisible/tetranucleotide_count/DNA_mechanics/PGLS/labelled_supertree_ottnames.tre")

## PREPROCESS DATA ----
matched_data <- data %>%
  filter(!is.na(Tip_Name), Tip_Name %in% phylo_tree$tip.label) %>%
  rename(species_name = Tip_Name)

phylo_tree <- compute.brlen(phylo_tree, method = "Grafen")
phylo_tree <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, matched_data$species_name))

## ---- Define which dinuc OE columns to use (collapsed for genome + intergenic) ----

collapsed_regions <- c("genome", "intergenic")

# Original 16 (non-collapsed) dinucleotide O/E columns (for other regions)
dinuc_oe_16 <- grep("_OE$", colnames(matched_data), value = TRUE)

# Remove mono-nucleotide OE columns if present
dinuc_oe_16 <- setdiff(dinuc_oe_16, c("A_OE", "T_OE", "G_OE", "C_OE"))

# Remove collapsed OE columns from the 16-set:
dinuc_oe_16 <- dinuc_oe_16[!grepl("_and_", dinuc_oe_16)]
dinuc_oe_16 <- dinuc_oe_16[!grepl("_OE_rc$", dinuc_oe_16)]

# Collapsed set (10 metrics): 6 paired collapsed + 4 palindromes
dinuc_oe_collapsed_10 <- c(
  grep("_and_.*_OE$", colnames(matched_data), value = TRUE),
  "AT_OE", "TA_OE", "CG_OE", "GC_OE"
) %>% unique()

dinuc_oe_collapsed_10 <- dinuc_oe_collapsed_10[dinuc_oe_collapsed_10 %in% colnames(matched_data)]

# Optional: stable order
dinuc_oe_16 <- sort(dinuc_oe_16)
dinuc_oe_collapsed_10 <- sort(dinuc_oe_collapsed_10)


## RUN PGLS WITH INTERACTION TERMS ----
pgls_results_all_regions <- list()

for (region in regions) {
  cat("Running PGLS with interaction terms for region:", region, "\n")

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
      ln_Genome_Size_c = as.numeric(scale(ln_Genome_Size, center = TRUE, scale = FALSE)),
      ln_lifespan_c    = as.numeric(scale(ln_lifespan, center = TRUE, scale = FALSE))
    )

  ## Choose motif set depending on region
  dinuc_columns_to_use <- if (region %in% collapsed_regions) dinuc_oe_collapsed_10 else dinuc_oe_16

  if (length(dinuc_columns_to_use) == 0) {
    message("No dinucleotide OE columns found for region: ", region)
    next
  }

  for (dinuc in dinuc_columns_to_use) {
    cat("  Running PGLS for", dinuc, "in region", region, "\n")

    region_df <- region_data %>%
      transmute(
        species_name     = species_name,
        ln_lifespan_c    = as.numeric(ln_lifespan_c),
        ln_Genome_Size_c = as.numeric(ln_Genome_Size_c),
        oe_raw           = suppressWarnings(as.numeric(.data[[dinuc]]))
      ) %>%
      # Must be > 0 to take ln
      filter(!is.na(oe_raw), oe_raw > 0) %>%
      mutate(
        ln_oe = log(oe_raw)
      ) %>%
      drop_na(ln_lifespan_c, ln_Genome_Size_c, ln_oe) %>%
      filter(is.finite(ln_lifespan_c), is.finite(ln_Genome_Size_c), is.finite(ln_oe))

    if (nrow(region_df) < 3) {
      message("    Skipping (complete cases < 3): ", dinuc, " in ", region, " (n=", nrow(region_df), ")")
      next
    }

    ## Model: ln(O/E) as response; centered ln(lifespan) and centered ln(genome size)
    formula_str <- "ln_oe ~ ln_lifespan_c * ln_Genome_Size_c"

    tryCatch({
      pgls_model <- gls(
        as.formula(formula_str),
        data = region_df,
        correlation = corPagel(1, phy = phylo_tree, form = ~species_name),
        method = "ML"
      )

      model_summary <- summary(pgls_model)$tTable

      ## Extract lambda
      cor_int <- tryCatch(intervals(pgls_model)$corStruct, error = function(e) NULL)
      lambda_value <- tryCatch(
        as.numeric(coef(pgls_model$modelStruct$corStruct, unconstrained = FALSE)),
        error = function(e) NA_real_
      )

      ## Term-name safe extraction (robust to ordering)
      get_term <- function(term, col) {
        if (term %in% rownames(model_summary)) model_summary[term, col] else NA_real_
      }

      term_lifespan <- "ln_lifespan_c"
      term_genome   <- "ln_Genome_Size_c"
      term_int_1    <- "ln_lifespan_c:ln_Genome_Size_c"
      term_int_2    <- "ln_Genome_Size_c:ln_lifespan_c"
      term_int      <- if (term_int_1 %in% rownames(model_summary)) term_int_1 else term_int_2

      pgls_results_all_regions[[paste(region, dinuc, sep = "_")]] <- data.frame(
        Region = region,
        Dinucleotide = dinuc,

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

        Coefficient_Interaction = get_term(term_int, "Value"),
        Std.Error_Interaction = get_term(term_int, "Std.Error"),
        t.value_Interaction = get_term(term_int, "t-value"),
        p.value_Interaction = get_term(term_int, "p-value"),

        Lambda = lambda_value,
        N_species = length(unique(region_df$species_name)),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      message("Error for ", dinuc, " in region ", region, ": ", e$message)
    })
  }
}

## COMBINE RESULTS + FDR ----
if (length(pgls_results_all_regions) > 0) {
  pgls_results <- bind_rows(pgls_results_all_regions) %>%
    mutate(
      p.adjusted.fdr_Intercept            = p.adjust(p.value_Intercept, method = "BH"),
      p.adjusted.fdr_ln_lifespan_c        = p.adjust(p.value_ln_lifespan_c, method = "BH"),
      p.adjusted.fdr_ln_Genome_Size_c     = p.adjust(p.value_ln_Genome_Size_c, method = "BH"),
      p.adjusted.fdr_Interaction          = p.adjust(p.value_Interaction, method = "BH")
    )

  out_file <- "pgls_results_8_regions_genome10_intergenic10_others16_dinuc_ln_OE_LS_CENTERED_G_CENTERED_interaction_OE_as_response.csv"
  write.csv(pgls_results, out_file, row.names = FALSE)
  message("Saved: ", out_file)
} else {
  message("No PGLS results generated.")
}

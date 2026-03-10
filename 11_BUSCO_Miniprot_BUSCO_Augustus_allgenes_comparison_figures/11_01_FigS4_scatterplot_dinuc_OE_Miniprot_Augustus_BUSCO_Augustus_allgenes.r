library(ggplot2)
library(dplyr)

# ============================================================
# INPUT FILES (new model form: ln(O/E) ~ ln(lifespan)c * ln(genome)c)
# ============================================================
file_miniprot_busco <- "pgls_results_8_regions_genome10_intergenic10_others16_dinuc_ln_OE_LS_CENTERED_G_CENTERED_interaction_OE_as_response_Miniprot_BUSCO.csv"
file_augustus_busco <- "pgls_results_8_regions_genome10_intergenic10_others16_dinuc_ln_OE_LS_CENTERED_G_CENTERED_interaction_OE_as_response_Augustus_BUSCO.csv"
file_augustus_all   <- "pgls_results_8_regions_genome10_intergenic10_others16_dinuc_ln_OE_LS_CENTERED_G_CENTERED_interaction_OE_as_response_Augustus_all_genes.csv"

# ============================================================
# READ
# ============================================================
miniprot_busco <- read.csv(file_miniprot_busco, check.names = FALSE)
augustus_busco <- read.csv(file_augustus_busco, check.names = FALSE)
augustus_all   <- read.csv(file_augustus_all,   check.names = FALSE)

# ============================================================
# BASIC CHECKS
# ============================================================
req <- c("Region", "Dinucleotide", "Coefficient_ln_lifespan_c", "Coefficient_ln_Genome_Size_c")
for (nm in c("miniprot_busco","augustus_busco","augustus_all")) {
  x <- get(nm)
  miss <- setdiff(req, colnames(x))
  if (length(miss) > 0) stop("Missing columns in ", nm, ": ", paste(miss, collapse = ", "))
}

# ============================================================
# Helper: make facet labels prettier (downstream1 -> Downstream 1)
# ============================================================
pretty_region <- function(x) {
  x <- gsub("([a-zA-Z]+)([0-9]+)", "\\1 \\2", x)
  tools::toTitleCase(x)
}

# ============================================================
# Helper: plotmath label for Sum|beta| (no unicode subscripts)
# beta_idx should be 1 or 2
# ============================================================
make_sum_label_plotmath <- function(sum_x, sum_y, label_x, label_y, beta_idx) {
  paste0(
    "atop(",
    "paste('Sum |', beta[", beta_idx, "], '| ", label_x, ": ', ", round(sum_x, 3), "),",
    "paste('Sum |', beta[", beta_idx, "], '| ", label_y, ": ', ", round(sum_y, 3), ")",
    ")"
  )
}

# ============================================================
# Core: build merged table for one coefficient column
#   coef_col = "Coefficient_ln_lifespan_c" or "Coefficient_ln_Genome_Size_c"
# ============================================================
build_merged <- function(coef_col) {

  m1 <- miniprot_busco %>%
    select(Region, Dinucleotide, all_of(coef_col)) %>%
    rename(beta_miniprot_busco = !!coef_col)

  m2 <- augustus_busco %>%
    select(Region, Dinucleotide, all_of(coef_col)) %>%
    rename(beta_augustus_busco = !!coef_col)

  m3 <- augustus_all %>%
    select(Region, Dinucleotide, all_of(coef_col)) %>%
    rename(beta_augustus_all = !!coef_col)

  merged <- m1 %>%
    inner_join(m2, by = c("Region", "Dinucleotide")) %>%
    inner_join(m3, by = c("Region", "Dinucleotide")) %>%
    filter(!Region %in% c("genome", "intergenic"))   # remove genome and intergenic regions

  merged
}

# ============================================================
# Helper: compute per-region sums for two columns and build plotmath labels
# ============================================================
make_beta_sums <- function(dat, col_x, col_y, label_x, label_y, beta_idx) {
  dat %>%
    group_by(Region) %>%
    summarise(
      sum_abs_x = sum(abs(.data[[col_x]]), na.rm = TRUE),
      sum_abs_y = sum(abs(.data[[col_y]]), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      label_plotmath = make_sum_label_plotmath(sum_abs_x, sum_abs_y, label_x, label_y, beta_idx)
    )
}

# ============================================================
# Helper: global symmetric axis limits for a given merged dataset
# ============================================================
get_axis_limits <- function(dat, cols, expand = 1.05) {
  vals <- unlist(dat[, cols], use.names = FALSE)
  lim <- max(abs(vals), na.rm = TRUE) * expand
  if (!is.finite(lim) || lim == 0) lim <- 1e-3
  c(-lim, lim)
}

# ============================================================
# Helper: draw a faceted scatterplot + save
# ============================================================
plot_faceted <- function(dat, xcol, ycol, sums_df, axis_limits,
                         x_lab_expr, y_lab_expr, out_png) {

  p <- ggplot(dat, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_point(size = 2.2, alpha = 0.70, color = "blue") +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    facet_wrap(
      ~ Region,
      scales = "fixed",
      labeller = labeller(Region = pretty_region)
    ) +
    geom_text(
      data = sums_df,
      aes(x = -Inf, y = Inf, label = label_plotmath),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1.1,
      size = 3.8,
      color = "red",
      parse = TRUE
    ) +
    scale_x_continuous(limits = axis_limits) +
    scale_y_continuous(limits = axis_limits) +
    coord_equal() +
    theme_minimal(base_size = 14) +
    labs(x = x_lab_expr, y = y_lab_expr)

  png(out_png, width = 1600, height = 1200, res = 150)
  print(p)
  dev.off()

  message("Saved: ", out_png)
}

# ============================================================
# SPECS for β1 and β2
# ============================================================
beta_specs <- list(
  list(
    coef_col = "Coefficient_ln_lifespan_c",
    beta_idx = 1,
    beta_tag = "beta1_lnLifespan",
    # axis label parts
    coef_text = "ln(lifespan)c"
  ),
  list(
    coef_col = "Coefficient_ln_Genome_Size_c",
    beta_idx = 2,
    beta_tag = "beta2_lnGenomeSize",
    coef_text = "ln(genome size)c"
  )
)

# ============================================================
# MAIN LOOP: make 4 plots (2 betas × 2 comparisons)
# ============================================================
for (sp in beta_specs) {

  merged <- build_merged(sp$coef_col)

  # shared limits across the three series, within this beta type
  axis_limits <- get_axis_limits(
    merged,
    cols = c("beta_miniprot_busco", "beta_augustus_busco", "beta_augustus_all"),
    expand = 1.05
  )

  # ----------------------------
  # Plot 1: Augustus all vs Augustus BUSCO
  # ----------------------------
  sums_all_vs_busco <- make_beta_sums(
    merged,
    col_x = "beta_augustus_all",
    col_y = "beta_augustus_busco",
    label_x = "Augustus all genes",
    label_y = "Augustus BUSCOs",
    beta_idx = sp$beta_idx
  )

  out1 <- paste0("scatterplots_by_region_", sp$beta_tag, "_Augustus_all_vs_BUSCOs.png")

  plot_faceted(
    dat = merged,
    xcol = "beta_augustus_all",
    ycol = "beta_augustus_busco",
    sums_df = sums_all_vs_busco,
    axis_limits = axis_limits,
    x_lab_expr = bquote(paste("PGLS coefficient (", beta[.(sp$beta_idx)], ") — Augustus all genes")),
    y_lab_expr = bquote(paste("PGLS coefficient (", beta[.(sp$beta_idx)], ") — Augustus BUSCOs")),
    out_png = out1
  )

  # ----------------------------
  # Plot 2: Miniprot BUSCO vs Augustus BUSCO
  # ----------------------------
  sums_miniprot_vs_aug <- make_beta_sums(
    merged,
    col_x = "beta_miniprot_busco",
    col_y = "beta_augustus_busco",
    label_x = "Miniprot BUSCOs",
    label_y = "Augustus BUSCOs",
    beta_idx = sp$beta_idx
  )

  out2 <- paste0("scatterplots_by_region_", sp$beta_tag, "_Miniprot_vs_Augustus_BUSCOs.png")

  plot_faceted(
    dat = merged,
    xcol = "beta_miniprot_busco",
    ycol = "beta_augustus_busco",
    sums_df = sums_miniprot_vs_aug,
    axis_limits = axis_limits,
    x_lab_expr = bquote(paste("PGLS coefficient (", beta[.(sp$beta_idx)], ") — Miniprot BUSCOs")),
    y_lab_expr = bquote(paste("PGLS coefficient (", beta[.(sp$beta_idx)], ") — Augustus BUSCOs")),
    out_png = out2
  )
}

message("Done.")

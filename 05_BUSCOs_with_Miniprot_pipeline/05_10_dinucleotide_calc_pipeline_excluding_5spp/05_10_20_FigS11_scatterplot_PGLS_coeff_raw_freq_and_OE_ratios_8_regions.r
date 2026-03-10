library(ggplot2)
library(dplyr)

# ============================================================
# INPUT FILES
# ============================================================
file_lnOE <- "pgls_results_8_regions_genome10_intergenic10_others16_dinuc_ln_OE_LS_CENTERED_G_CENTERED_interaction_OE_as_response.csv"
file_freq <- "pgls_results_8_regions_genome10_intergenic10_others16_dinuc_RATE_LS_CENTERED_G_CENTERED_interaction_raw_freq_as_response.csv"

# ============================================================
# REGIONS TO PLOT (8 panels; keep this order)
# ============================================================
regions_to_keep <- c("genome", "intergenic", "upstream1", "upstream2", "exon", "intron", "downstream1", "downstream2")

# ============================================================
# READ
# ============================================================
dat_lnOE <- read.csv(file_lnOE, check.names = FALSE)
dat_freq <- read.csv(file_freq, check.names = FALSE)

# ============================================================
# FACET LABELS (capitalize + add space before digits)
# ============================================================
facet_title <- function(x) {
  x <- as.character(x)
  x <- gsub("([a-zA-Z]+)([0-9]+)", "\\1 \\2", x)  # upstream1 -> upstream 1
  tools::toTitleCase(x)
}

# ============================================================
# PREP: shared dinucleotide key
#   - freq Dinucleotide: "AA_Rate" or "AA_and_TT_Rate" -> "AA" or "AA_and_TT"
#   - lnOE Dinucleotide: "AA_OE"   or "AA_and_TT_OE"   -> "AA" or "AA_and_TT"
# NOTE: We intentionally keep the collapsed labels as collapsed (AA_and_TT, etc.),
#       so genome/intergenic join correctly on the collapsed metrics.
# ============================================================
dat_lnOE2_ls <- dat_lnOE %>%
  select(Region, Dinucleotide, Coefficient_ln_lifespan_c) %>%
  mutate(
    Region = as.character(Region),
    Dinucleotide = as.character(Dinucleotide),
    Dinuc = sub("_OE$", "", Dinucleotide)
  ) %>%
  rename(beta_lnOE_ls = Coefficient_ln_lifespan_c)

dat_lnOE2_gs <- dat_lnOE %>%
  select(Region, Dinucleotide, Coefficient_ln_Genome_Size_c) %>%
  mutate(
    Region = as.character(Region),
    Dinucleotide = as.character(Dinucleotide),
    Dinuc = sub("_OE$", "", Dinucleotide)
  ) %>%
  rename(beta_lnOE_gs = Coefficient_ln_Genome_Size_c)

dat_freq2_ls <- dat_freq %>%
  select(Region, Dinucleotide, Coefficient_ln_lifespan_c) %>%
  mutate(
    Region = as.character(Region),
    Dinucleotide = as.character(Dinucleotide),
    Dinuc = sub("_Rate$", "", Dinucleotide)
  ) %>%
  rename(beta_freq_ls = Coefficient_ln_lifespan_c)

dat_freq2_gs <- dat_freq %>%
  select(Region, Dinucleotide, Coefficient_ln_Genome_Size_c) %>%
  mutate(
    Region = as.character(Region),
    Dinucleotide = as.character(Dinucleotide),
    Dinuc = sub("_Rate$", "", Dinucleotide)
  ) %>%
  rename(beta_freq_gs = Coefficient_ln_Genome_Size_c)

# ============================================================
# MERGE: lifespan coefficients (β1)
# ============================================================
dat_m_ls <- dat_freq2_ls %>%
  inner_join(dat_lnOE2_ls, by = c("Region", "Dinuc")) %>%
  filter(
    Region %in% regions_to_keep,
    !is.na(beta_freq_ls),
    !is.na(beta_lnOE_ls)
  ) %>%
  mutate(Region = factor(Region, levels = regions_to_keep))

if (nrow(dat_m_ls) == 0) {
  stop("dat_m_ls has 0 rows after join+filter. Check Dinucleotide naming and columns.")
}

# ============================================================
# MERGE: genome-size coefficients (β2)
# ============================================================
dat_m_gs <- dat_freq2_gs %>%
  inner_join(dat_lnOE2_gs, by = c("Region", "Dinuc")) %>%
  filter(
    Region %in% regions_to_keep,
    !is.na(beta_freq_gs),
    !is.na(beta_lnOE_gs)
  ) %>%
  mutate(Region = factor(Region, levels = regions_to_keep))

if (nrow(dat_m_gs) == 0) {
  stop("dat_m_gs has 0 rows after join+filter. Check Dinucleotide naming and columns.")
}

# ============================================================
# PER-REGION STATS: Pearson r + Bonferroni-adjusted p (across 8 panels)
# ============================================================
stats_by_region <- function(dat, xcol, ycol) {
  dat %>%
    group_by(Region) %>%
    summarise(
      r = suppressWarnings(cor(.data[[xcol]], .data[[ycol]], method = "pearson", use = "complete.obs")),
      p = suppressWarnings(cor.test(.data[[xcol]], .data[[ycol]], method = "pearson")$p.value),
      n = sum(is.finite(.data[[xcol]]) & is.finite(.data[[ycol]])),
      .groups = "drop"
    ) %>%
    mutate(
      adj_p = p.adjust(p, method = "BH"),
      label = paste0(
        "r = ", formatC(r, format = "f", digits = 2), "\n",
        "adj. p = ", formatC(adj_p, format = "e", digits = 2), "\n",
        "n = ", n
      )
    )
}

# ============================================================
# PLOT FUNCTION
# ============================================================
plot_compare <- function(dat, xcol, ycol, stats_df, out_png, xlab_expr, ylab_expr) {

  # Common axis limits across ALL panels for this plot
  x_max <- max(abs(dat[[xcol]]), na.rm = TRUE)
  y_max <- max(abs(dat[[ycol]]), na.rm = TRUE)
  x_limits <- c(-x_max, x_max)
  y_limits <- c(-y_max, y_max)

  png(out_png, width = 2000, height = 1400, res = 150)

  p <- ggplot(dat, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_point(size = 2.0, color = "black", alpha = 0.70) +
    geom_smooth(method = "lm", se = FALSE,
                color = adjustcolor("black", alpha.f = 0.35),
                linewidth = 1) +
    facet_wrap(
      ~ Region,
      scales = "fixed",
      labeller = labeller(Region = facet_title)
    ) +
    geom_text(
      data = stats_df,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.05, vjust = 1.1,
      size = 3.4,
      color = "black",
      lineheight = 1.05
    ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits) +
    theme_minimal(base_size = 14) +
    theme(
      text = element_text(family = "Arial", colour = "black"),
      strip.text = element_text(size = 13, face = "bold", colour = "black"),
      axis.title = element_text(size = 14, colour = "black"),
      axis.text  = element_text(size = 12, colour = "black")
    ) +
    labs(x = xlab_expr, y = ylab_expr)

  print(p)
  dev.off()

  cat("Saved:", out_png, "\n")
}

# ============================================================
# 1) Lifespan coefficients comparison (β1)
#   x: raw frequency model β1 (lifespan term)
#   y: ln(O/E) model        β1 (lifespan term)
# ============================================================
stats_ls <- stats_by_region(dat_m_ls, "beta_freq_ls", "beta_lnOE_ls")

plot_compare(
  dat = dat_m_ls,
  xcol = "beta_freq_ls",
  ycol = "beta_lnOE_ls",
  stats_df = stats_ls,
  out_png = "scatterplots_by_region_beta1_rawfreq_vs_lnOE_8regions.png",
  xlab_expr = expression(paste("PGLS coefficient (", beta[1], ") — raw dinucleotide frequency")),
  ylab_expr = expression(paste("PGLS coefficient (", beta[1], ") — ln(O/E)"))
)

# ============================================================
# 2) Genome size coefficients comparison (β2)
#   x: raw frequency model β2 (genome size term)
#   y: ln(O/E) model        β2 (genome size term)
# ============================================================
stats_gs <- stats_by_region(dat_m_gs, "beta_freq_gs", "beta_lnOE_gs")

plot_compare(
  dat = dat_m_gs,
  xcol = "beta_freq_gs",
  ycol = "beta_lnOE_gs",
  stats_df = stats_gs,
  out_png = "scatterplots_by_region_beta2_rawfreq_vs_lnOE_8regions.png",
  xlab_expr = expression(paste("PGLS coefficient (", beta[2], ") — raw dinucleotide frequency")),
  ylab_expr = expression(paste("PGLS coefficient (", beta[2], ") — ln(O/E)"))
)

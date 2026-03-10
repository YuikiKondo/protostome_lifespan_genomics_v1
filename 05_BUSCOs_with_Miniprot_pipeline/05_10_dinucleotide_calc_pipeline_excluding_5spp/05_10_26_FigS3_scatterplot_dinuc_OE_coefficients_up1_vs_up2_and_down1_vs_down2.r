# Compare PGLS coefficients between:
#  - upstream1 vs upstream2
#  - downstream1 vs downstream2
# For BOTH coefficient types:
#   (1) Lifespan coefficient  β₁  : Coefficient_ln_lifespan_c
#   (2) Genome size coefficient β₂: Coefficient_ln_Genome_Size_c
#
# Features:
#  - SAME x/y range for upstream and downstream within each coefficient type (tight symmetric limit)
#  - Larger text everywhere
#  - Points + regression line
#  - Annotation: regression equation + Sum |β| for each region
#  - Transparent annotation background (doesn't hide points)
#  - Axis title wording switches appropriately:
#       β₁ -> "Lifespan association"
#       β₂ -> "Genome size association"
#  - No diagonal line, no point labels

library(tidyverse)

# ---- INPUT ----
input_file <- "pgls_results_8_regions_genome10_intergenic10_others16_dinuc_ln_OE_LS_CENTERED_G_CENTERED_interaction_OE_as_response.csv"
df <- read.csv(input_file, check.names = FALSE)

# ---- CHECK FORMAT ----
required_cols <- c("Region", "Dinucleotide", "Coefficient_ln_lifespan_c", "Coefficient_ln_Genome_Size_c")
missing <- setdiff(required_cols, colnames(df))
if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

regions_of_interest <- c("upstream1", "upstream2", "downstream1", "downstream2")

# ---- helper: round limit to a "nice" small step (keeps small ranges tight) ----
nice_lim <- function(x, margin = 0.08) {
  x <- as.numeric(x)
  if (!is.finite(x) || x <= 0) return(0.05)

  x <- x * (1 + margin)

  step <- if (x < 0.05) 0.005 else if (x < 0.2) 0.01 else if (x < 1) 0.05 else 0.1
  ceiling(x / step) * step
}

# ---- PLOTTING FUNCTION (generic for any coefficient column) ----
plot_region_pair <- function(df, region_x, region_y,
                             coef_col,
                             out_png,
                             beta_symbol = "β",
                             assoc_label = "Lifespan association",
                             lim) {

  dat <- df %>%
    filter(Region %in% c(region_x, region_y)) %>%
    select(Dinucleotide, Region, all_of(coef_col)) %>%
    mutate(
      Region = as.character(Region),
      coef = suppressWarnings(as.numeric(.data[[coef_col]]))
    ) %>%
    select(Dinucleotide, Region, coef) %>%
    pivot_wider(names_from = Region, values_from = coef) %>%
    filter(is.finite(.data[[region_x]]), is.finite(.data[[region_y]]))

  if (nrow(dat) == 0) {
    stop("No finite pairs found for: ", region_x, " vs ", region_y, " using ", coef_col)
  }

  # Sum |β| for each region (within the dinucleotides present)
  sum_abs_x <- sum(abs(dat[[region_x]]), na.rm = TRUE)
  sum_abs_y <- sum(abs(dat[[region_y]]), na.rm = TRUE)

  # Regression (region_y ~ region_x)
  lm_fit <- lm(dat[[region_y]] ~ dat[[region_x]])
  b0 <- unname(coef(lm_fit)[1])
  b1 <- unname(coef(lm_fit)[2])

  # Annotation text (transparent label)
  ann_text <- paste(
    sprintf("y = %.3fx %+.3f", b1, b0),
    sprintf("Sum |%s| (%s): %.3f", beta_symbol, region_x, sum_abs_x),
    sprintf("Sum |%s| (%s): %.3f", beta_symbol, region_y, sum_abs_y),
    sep = "\n"
  )

  # Axis labels (switch wording for β₁ vs β₂ via assoc_label)
  x_label <- paste0(
    assoc_label, " of ln(O/E)\n(",
    beta_symbol, ") in ", tools::toTitleCase(region_x)
  )
  y_label <- paste0(
    assoc_label, " of ln(O/E)\n(",
    beta_symbol, ") in ", tools::toTitleCase(region_y)
  )

  p <- ggplot(dat, aes(x = .data[[region_x]], y = .data[[region_y]])) +
    geom_point(size = 4) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.4) +

    annotate(
      "label",
      x = -Inf, y = Inf,
      label = ann_text,
      hjust = -0.05, vjust = 1.1,
      size = 7,
      label.size = 0,
      fill = NA
    ) +

    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim), clip = "off") +
    labs(title = NULL, x = x_label, y = y_label) +
    theme_minimal(base_size = 18) +
    theme(
      axis.title = element_text(size = 24),
      axis.text  = element_text(size = 24),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.3),
      plot.margin = margin(10, 18, 10, 28)
    )

  ggsave(out_png, p, width = 9, height = 8, dpi = 300)
  p
}

# ============================================================
# Plot specs for β₁ and β₂
# ============================================================
coef_specs <- list(
  list(
    coef_col = "Coefficient_ln_lifespan_c",
    beta_symbol = "β\u2081",            # β₁
    assoc_label = "Lifespan association",
    tag = "beta1_lnLifespan"
  ),
  list(
    coef_col = "Coefficient_ln_Genome_Size_c",
    beta_symbol = "β\u2082",            # β₂
    assoc_label = "Genome size association",
    tag = "beta2_lnGenomeSize"
  )
)

# ============================================================
# Make 4 plots (2 region-pairs × 2 coefficient types)
# ============================================================
for (spec in coef_specs) {

  # tight global limit for this coefficient type (shared across upstream+downstream plots)
  global_max <- df %>%
    filter(Region %in% regions_of_interest) %>%
    mutate(val = suppressWarnings(as.numeric(.data[[spec$coef_col]]))) %>%
    pull(val) %>%
    abs() %>%
    max(na.rm = TRUE)

  lim <- nice_lim(global_max, margin = 0.08)

  # ---- Upstream: upstream1 vs upstream2 ----
  p_up <- plot_region_pair(
    df,
    region_x = "upstream1",
    region_y = "upstream2",
    coef_col = spec$coef_col,
    out_png  = paste0("upstream2_vs_upstream1_", spec$tag, "_tight_range.png"),
    beta_symbol = spec$beta_symbol,
    assoc_label = spec$assoc_label,
    lim = lim
  )
  print(p_up)

  # ---- Downstream: downstream1 vs downstream2 ----
  p_down <- plot_region_pair(
    df,
    region_x = "downstream1",
    region_y = "downstream2",
    coef_col = spec$coef_col,
    out_png  = paste0("downstream2_vs_downstream1_", spec$tag, "_tight_range.png"),
    beta_symbol = spec$beta_symbol,
    assoc_label = spec$assoc_label,
    lim = lim
  )
  print(p_down)

  message("Saved plots for ", spec$coef_col, " with global symmetric limit ±", lim)
}

message("Done.")

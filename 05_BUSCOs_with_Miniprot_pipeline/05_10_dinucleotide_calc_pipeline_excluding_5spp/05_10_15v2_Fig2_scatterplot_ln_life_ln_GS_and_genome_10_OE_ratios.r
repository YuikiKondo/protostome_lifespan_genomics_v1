## ============================================================
## Fig2A & Fig2B: 10-panel scatterplots (genome rows only)
##   Fig2A: x = ln(lifespan [days])
##   Fig2B: x = ln(genome size [bp])
## y = dinucleotide O/E (10 collapsed metrics, shared y-scale)
## Add: OLS equation + R^2 (NO p-value)
## Output: 2 PNGs
## ============================================================

## LOAD LIBRARIES ----
library(tidyverse)
library(ggh4x)
library(grid)  # unit()

## LOAD DATA ----
input_file <- "di_obs_exp_ratio_rev_comp_collapsed_only_with_lifespan_tip_names_and_genome_size_Github_check.csv"
data <- read.csv(input_file, check.names = FALSE)

## BASIC PREP ----
data$Average_lifespan_days <- suppressWarnings(as.numeric(data$Average_lifespan_days))
data$ln_lifespan <- ifelse(
  !is.na(data$Average_lifespan_days) & data$Average_lifespan_days > 0,
  log(data$Average_lifespan_days),
  NA_real_
)

data$Genome_size_bp <- suppressWarnings(as.numeric(data$Genome_size_bp))
data$ln_genome_size <- ifelse(
  !is.na(data$Genome_size_bp) & data$Genome_size_bp > 0,
  log(data$Genome_size_bp),
  NA_real_
)

## DINUCLEOTIDE COLUMNS (fixed order) ----
motif_cols <- c(
  "AT_OE", "TA_OE", "GC_OE", "CG_OE",
  "CA_and_TG_OE", "AC_and_GT_OE", "AG_and_CT_OE", "GA_and_TC_OE",
  "AA_and_TT_OE", "CC_and_GG_OE"
)

motif_labels <- c(
  AT_OE = "AT",
  TA_OE = "TA",
  GC_OE = "GC",
  CG_OE = "CG",
  CA_and_TG_OE = "CA + TG",
  AC_and_GT_OE = "AC + GT",
  AG_and_CT_OE = "AG + CT",
  GA_and_TC_OE = "GA + TC",
  AA_and_TT_OE = "AA + TT",
  CC_and_GG_OE = "CC + GG"
)

## FILTER: genome region only + coerce motif cols to numeric ----
d <- data %>%
  filter(Region == "genome") %>%
  mutate(across(all_of(motif_cols), ~ suppressWarnings(as.numeric(.x))))

missing_cols <- setdiff(motif_cols, colnames(d))
if (length(missing_cols) > 0) {
  stop("Missing motif columns: ", paste(missing_cols, collapse = ", "))
}

## LONG FORMAT (keep BOTH x variables) ----
plot_df <- d %>%
  select(ln_lifespan, ln_genome_size, all_of(motif_cols)) %>%
  pivot_longer(
    cols = all_of(motif_cols),
    names_to = "Motif",
    values_to = "OE"
  ) %>%
  drop_na(OE) %>%
  mutate(Motif = factor(Motif, levels = motif_cols))

## GLOBAL Y-LIMITS (shared across panels; use all available points) ----
y_limits <- range(plot_df$OE, na.rm = TRUE)

## ------------------------------------------------------------
## Helper: build per-facet annotation df (equation + R^2 only)
## ------------------------------------------------------------
make_lm_labels <- function(df_long, x_col, y_col = "OE") {

  df_use <- df_long %>% filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]))

  stats_df <- df_use %>%
    group_by(Motif) %>%
    group_modify(~{
      fit <- lm(reformulate(x_col, y_col), data = .x)
      sm  <- summary(fit)

      a <- unname(coef(fit)[1])
      b <- unname(coef(fit)[2])
      r2 <- unname(sm$r.squared)

      # Keep formatting compact & consistent
      label <- sprintf("y = %.3f %s %.3f x\nR² = %.3f",
                       a,
                       ifelse(b >= 0, "+", "−"),
                       abs(b),
                       r2)

      tibble(label = label)
    }) %>%
    ungroup()

  # Position text consistently in each facet using panel-wise ranges
  pos_df <- df_use %>%
    group_by(Motif) %>%
    summarise(
      x_min = min(.data[[x_col]], na.rm = TRUE),
      x_max = max(.data[[x_col]], na.rm = TRUE),
      y_min = min(.data[[y_col]], na.rm = TRUE),
      y_max = max(.data[[y_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      x = x_min + 0.05 * (x_max - x_min),
      y = y_max - 0.05 * (y_max - y_min)
    ) %>%
    select(Motif, x, y)

  left_join(stats_df, pos_df, by = "Motif")
}

## ------------------------------------------------------------
## Plot helper: make 10-panel plot for a chosen x variable
## ------------------------------------------------------------
make_panel_plot <- function(df_long, x_col, x_lab, out_png) {

  df_use <- df_long %>%
    filter(!is.na(.data[[x_col]]))

  ann_df <- make_lm_labels(df_use, x_col = x_col, y_col = "OE")

  p <- ggplot(df_use, aes(x = .data[[x_col]], y = OE)) +
    geom_point(size = 0.9, alpha = 0.45) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      colour = "black",
      linewidth = 1.0,
      alpha = 0.6
    ) +
    geom_text(
      data = ann_df,
      aes(label = label),
      x = -Inf,
      y = Inf,
      inherit.aes = FALSE,
      hjust = -0.1,   # slightly inside from left border
      vjust = 1.25,    # slightly below top border
      size = 4.2,
      family = "Arial",
      colour = "black",
      lineheight = 0.95
    ) +
    facet_wrap2(
      ~ Motif,
      ncol = 5,
      axes = "all",
      remove_labels = "all",
      labeller = labeller(Motif = motif_labels)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.10))) +
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.10))) +
    coord_cartesian(ylim = y_limits, clip = "off") +
    theme_classic(base_size = 18) +
    theme(
      text = element_text(family = "Arial", colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(size = 14, face = "bold", colour = "black"),

      axis.title.x = element_text(size = 20, margin = margin(t = 12), colour = "black"),
      axis.title.y = element_text(size = 20, margin = margin(r = 12), colour = "black"),

      axis.text = element_text(size = 16, colour = "black"),

      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks.length = unit(6, "pt"),

      axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0)
    ) +
    labs(
      x = x_lab,
      y = "Dinucleotide O/E ratio"
    )

  ggsave(out_png, p, width = 16, height = 7.5, dpi = 300)
  cat("Saved: ", out_png, "\n", sep = "")
}

## ------------------------------------------------------------
## Fig2A: ln(lifespan) vs O/E
## ------------------------------------------------------------
make_panel_plot(
  df_long = plot_df,
  x_col   = "ln_lifespan",
  x_lab   = "ln(Lifespan [days])",
  out_png = "Fig2A_genome_dinuc_OE_vs_lnLifespan_10panels.png"
)

## ------------------------------------------------------------
## Fig2B: ln(genome size) vs O/E
## ------------------------------------------------------------
make_panel_plot(
  df_long = plot_df,
  x_col   = "ln_genome_size",
  x_lab   = "ln(Genome size [bp])",
  out_png = "Fig2B_genome_dinuc_OE_vs_lnGenomeSize_10panels.png"
)
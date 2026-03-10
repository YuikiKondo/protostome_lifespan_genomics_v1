## ============================================================
## SupFig: Intergenic proportion vs ln(lifespan) and ln(genome size)
##
## GOAL
##   (A) Univariate PGLS: intergenic_prop ~ ln(lifespan)
##   (B) Univariate PGLS: intergenic_prop ~ ln(genome size)
##   (C) Multivariate PGLS (centered): intergenic_prop ~ ln(lifespan)c * ln(genome)c
##       + visualize interaction as marginal-effect lines (Q1/Q2/Q3 lifespan)
##
## KEY POINTS
##   - Use genome rows only (Region == "genome")
##   - Intergenic proportion = Intergenic_bp / Genome_size_bp
##   - PGLS via nlme::gls with corPagel (ML), Grafen branch lengths
##   - Consistent styling (Arial, black, transparent annotations)
##   - Panel C legend: forced to 2 rows + slightly taller device to avoid clipping
## ============================================================

library(tidyverse)
library(grid)       # unit()
library(ape)
library(nlme)
library(patchwork)

## ---- INPUT ----
input_file <- "filtered_species_with_tip_name_dinuc_Github_check_plus_intergenic.csv"
df_raw <- read.csv(input_file, check.names = FALSE)

## ---- CHECK REQUIRED COLUMNS ----
required_cols <- c("Region", "Tip_Name", "Average_lifespan_days", "Genome_size_bp",
                   "Intergenic_bp", "Total_bp_gene_associated_region_including_N")
missing <- setdiff(required_cols, colnames(df_raw))
if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

## ============================================================
## 1) PREP DATA (genome rows only + derived variables)
## ============================================================

df <- df_raw %>%
  mutate(
    Region = as.character(Region),
    species_name = as.character(Tip_Name),

    # coerce numeric safely
    Average_lifespan_days = suppressWarnings(as.numeric(Average_lifespan_days)),
    Genome_size_bp        = suppressWarnings(as.numeric(Genome_size_bp)),
    Intergenic_bp         = suppressWarnings(as.numeric(Intergenic_bp))
  ) %>%
  filter(Region == "genome") %>%
  mutate(
    # log transforms
    ln_lifespan = ifelse(!is.na(Average_lifespan_days) & Average_lifespan_days > 0,
                         log(Average_lifespan_days), NA_real_),
    ln_genome_size = ifelse(!is.na(Genome_size_bp) & Genome_size_bp > 0,
                            log(Genome_size_bp), NA_real_),

    # intergenic proportion
    intergenic_prop = ifelse(!is.na(Intergenic_bp) & !is.na(Genome_size_bp) & Genome_size_bp > 0,
                             Intergenic_bp / Genome_size_bp, NA_real_)
  )

if (nrow(df) == 0) stop("No genome rows found (Region == 'genome').")

if (any(df$intergenic_prop < 0 | df$intergenic_prop > 1, na.rm = TRUE)) {
  warning("Some intergenic_prop values outside [0,1]. Check Intergenic_bp and Genome_size_bp.")
}

## ============================================================
## 2) LOAD TREE + MATCH SPECIES
## ============================================================

phylo_tree <- read.tree("/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_exon_intron_separated/gene_region_separated_sequences_no_internal_stop_3divisible/tetranucleotide_count/DNA_mechanics/PGLS/labelled_supertree_ottnames.tre")
phylo_tree <- compute.brlen(phylo_tree, method = "Grafen")

# keep only species present in tree
df <- df %>%
  filter(!is.na(species_name), species_name %in% phylo_tree$tip.label)

if (nrow(df) < 3) stop("Too few genome species after matching to tree (n < 3).")

# prune tree to match df
phylo_tree <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, df$species_name))

## ============================================================
## 3) HELPER: Fit univariate PGLS y ~ x
##    Returns: data used, n, intercept, slope, p, lambda
## ============================================================

fit_pgls <- function(dat, x, y, tree) {

  d <- dat %>%
    transmute(
      species_name = species_name,
      x = suppressWarnings(as.numeric(.data[[x]])),
      y = suppressWarnings(as.numeric(.data[[y]]))
    ) %>%
    drop_na(species_name, x, y) %>%
    filter(is.finite(x), is.finite(y)) %>%
    filter(species_name %in% tree$tip.label)

  n <- nrow(d)
  if (n < 3) {
    return(list(df = d, n = n, intercept = NA_real_, slope = NA_real_, p = NA_real_, lambda = NA_real_))
  }

  tr <- drop.tip(tree, setdiff(tree$tip.label, d$species_name))

  fit <- gls(
    y ~ x,
    data = d,
    correlation = corPagel(1, phy = tr, form = ~species_name),
    method = "ML"
  )

  tt <- summary(fit)$tTable
  out <- list(
    df = d,
    n = n,
    intercept = unname(tt["(Intercept)", "Value"]),
    slope     = unname(tt["x", "Value"]),
    p         = unname(tt["x", "p-value"]),
    lambda    = tryCatch(as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE)),
                         error = function(e) NA_real_)
  )
  out
}

## ============================================================
## 4) HELPER: plot univariate panel with PGLS line + annotation
## ============================================================

plot_panel <- function(fit_obj, xlab, ylab, model_label = NULL,
                       x_fixed = NULL, y_fixed = c(0, 1.05)) {

  d <- fit_obj$df

  header <- if (!is.null(model_label)) paste0("Model: ", model_label, "\n") else ""
  ann <- paste0(
    header,
    "PGLS slope = ", formatC(fit_obj$slope, format = "f", digits = 3), "\n",
    " p = ", formatC(fit_obj$p, format = "e", digits = 2), "\n",
    "  λ = ", formatC(fit_obj$lambda, format = "f", digits = 2), "\n",
    "  n = ", fit_obj$n
  )

  # x limits with small padding
  if (is.null(x_fixed)) {
    x_rng <- range(d$x, na.rm = TRUE)
    x_pad <- diff(x_rng) * 0.05
    x_fixed <- c(x_rng[1] - x_pad, x_rng[2] + x_pad)
  }

  ggplot(d, aes(x = x, y = y)) +
    geom_point(size = 1.2, alpha = 0.50, color = "black") +
    geom_abline(intercept = fit_obj$intercept, slope = fit_obj$slope,
                colour = "black", linewidth = 1.2, alpha = 0.5) +

    # transparent annotation (does not hide points)
    geom_label(
      data = tibble(x = -Inf, y = Inf, label = ann),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = -0.05, vjust = 1.05,
      size = 3.8, color = "black",
      fill = NA, label.size = 0, label.r = unit(0, "pt"),
      lineheight = 1.05
    ) +

    coord_cartesian(xlim = x_fixed, ylim = y_fixed, clip = "off") +
    scale_y_continuous(
      breaks = c(0.00, 0.25, 0.50, 0.75, 1.00),
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +

    theme_classic(base_size = 16) +
    theme(
      text = element_text(family = "Arial", colour = "black"),
      axis.title.x = element_text(size = 18, margin = margin(t = 8)),
      axis.title.y = element_text(size = 18, margin = margin(r = 10)),
      axis.text = element_text(size = 14),
      axis.line = element_blank(),
      axis.ticks = element_line(linewidth = 0.4),
      axis.ticks.length = unit(5, "pt"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.3),
      plot.margin = margin(10, 18, 10, 28)
    ) +
    labs(x = xlab, y = ylab)
}

## ============================================================
## 5) PANELS A & B (univariate PGLS) + save combined figure
## ============================================================

# A: intergenic_prop ~ ln(lifespan)
fit_a <- fit_pgls(df, x = "ln_lifespan", y = "intergenic_prop", tree = phylo_tree)
pA <- plot_panel(
  fit_a,
  xlab = "ln(Lifespan [days])",
  ylab = "Intergenic proportion",
  model_label = "intergenic prop. ~ ln(lifespan)",
  y_fixed = c(0, 1.1)   # keep your preferred headroom
)

# B: intergenic_prop ~ ln(genome size)
fit_b <- fit_pgls(df, x = "ln_genome_size", y = "intergenic_prop", tree = phylo_tree)
pB <- plot_panel(
  fit_b,
  xlab = "ln(Genome size [bp])",
  ylab = "Intergenic proportion",
  model_label = "intergenic prop. ~ ln(genome size)",
  y_fixed = c(0, 1.1)
)

pAB <- pA | pB
ggsave("SupFig_intergenic_prop_vs_lnLifespan_lnGenome_PGLS.png",
       pAB, width = 12.5, height = 6.0, dpi = 300)
cat("Saved: SupFig_intergenic_prop_vs_lnLifespan_lnGenome_PGLS.png\n")

## ============================================================
## 6) MULTIVARIATE PGLS (centered) + export coefficient table
## ============================================================

cat("\n================ MULTIVARIATE PGLS (CENTERED) ================\n")
cat("Model: intergenic_prop ~ ln_lifespan_c * ln_genome_size_c\n")

d_multi <- df %>%
  transmute(
    species_name = species_name,
    intergenic_prop = suppressWarnings(as.numeric(intergenic_prop)),
    ln_lifespan = suppressWarnings(as.numeric(ln_lifespan)),
    ln_genome_size = suppressWarnings(as.numeric(ln_genome_size))
  ) %>%
  drop_na(species_name, intergenic_prop, ln_lifespan, ln_genome_size) %>%
  filter(is.finite(intergenic_prop), is.finite(ln_lifespan), is.finite(ln_genome_size)) %>%
  filter(species_name %in% phylo_tree$tip.label)

if (nrow(d_multi) < 3) stop("Too few rows for multivariate PGLS after filtering (n < 3).")

# center predictors (global centering)
d_multi <- d_multi %>%
  mutate(
    ln_lifespan_c = as.numeric(scale(ln_lifespan, center = TRUE, scale = FALSE)),
    ln_genome_size_c = as.numeric(scale(ln_genome_size, center = TRUE, scale = FALSE))
  )

tr_multi <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, d_multi$species_name))

fit_multi <- gls(
  intergenic_prop ~ ln_lifespan_c * ln_genome_size_c,
  data = d_multi,
  correlation = corPagel(1, phy = tr_multi, form = ~species_name),
  method = "ML"
)

tt <- summary(fit_multi)$tTable
lambda_multi <- tryCatch(as.numeric(coef(fit_multi$modelStruct$corStruct, unconstrained = FALSE)),
                         error = function(e) NA_real_)

cat("n =", nrow(d_multi), "\n")
cat("lambda =", formatC(lambda_multi, format = "f", digits = 3), "\n\n")

coef_df <- data.frame(
  term = rownames(tt),
  estimate = tt[, "Value"],
  se = tt[, "Std.Error"],
  t = tt[, "t-value"],
  p = tt[, "p-value"],
  lambda = lambda_multi,
  n = nrow(d_multi),
  row.names = NULL,
  check.names = FALSE
)
print(coef_df)

write.csv(coef_df,
          "PGLS_multivariate_intergenic_prop_lnLifespan_lnGenome_CENTERED_interaction.csv",
          row.names = FALSE)

cat("Wrote: PGLS_multivariate_intergenic_prop_lnLifespan_lnGenome_CENTERED_interaction.csv\n")
cat("=============================================================\n")

## ============================================================
## 7) PANEL C: visualize interaction as marginal-effect lines
##    (intergenic_prop vs ln_genome_size_c at lifespan Q1/Q2/Q3)
## ============================================================

# Choose lifespan slices (centered ln(lifespan) quantiles)
L_vals <- as.numeric(quantile(d_multi$ln_lifespan_c, probs = c(0.25, 0.50, 0.75), na.rm = TRUE))
names(L_vals) <- c("Short lifespan (Q1)", "Median lifespan (Q2)", "Long lifespan (Q3)")

# Genome-size grid (centered)
G_seq <- seq(min(d_multi$ln_genome_size_c, na.rm = TRUE),
             max(d_multi$ln_genome_size_c, na.rm = TRUE),
             length.out = 200)

# Prediction grid for 3 lifespan slices
pred_grid <- bind_rows(lapply(names(L_vals), function(lbl) {
  tibble(
    ln_genome_size_c = G_seq,
    ln_lifespan_c = L_vals[[lbl]],
    group = lbl
  )
}))

pred_grid$pred <- as.numeric(predict(fit_multi, newdata = pred_grid))

# Order legend lines (as displayed)
pred_grid$group <- factor(pred_grid$group,
                          levels = c("Long lifespan (Q3)", "Median lifespan (Q2)", "Short lifespan (Q1)"))

# Annotation for Panel C (full model + key terms)
ttm <- summary(fit_multi)$tTable
annC <- paste0(
  "Model: intergenic prop. ~ ln(lifespan)c * ln(genome)c\n",
  "β_L = ",  formatC(ttm["ln_lifespan_c","Value"], format="f", digits=3),
  ", p = ",  formatC(ttm["ln_lifespan_c","p-value"], format="e", digits=2), "\n",
  "β_G = ",  formatC(ttm["ln_genome_size_c","Value"], format="f", digits=3),
  ", p = ",  formatC(ttm["ln_genome_size_c","p-value"], format="e", digits=2), "\n",
  "β_L×G = ",formatC(ttm["ln_lifespan_c:ln_genome_size_c","Value"], format="f", digits=3),
  ", p = ",  formatC(ttm["ln_lifespan_c:ln_genome_size_c","p-value"], format="e", digits=2), "\n",
  "λ = ",    formatC(lambda_multi, format="f", digits=2),
  ", n = ",  nrow(d_multi)
)

pC <- ggplot(d_multi, aes(x = ln_genome_size_c, y = intergenic_prop)) +
  geom_point(size = 1.2, alpha = 0.35, color = "black") +
  geom_line(data = pred_grid, aes(y = pred, linetype = group),
            color = "black", linewidth = 1.1) +

  geom_label(
    data = tibble(x = -Inf, y = Inf, label = annC),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.05,
    size = 3.6, color = "black",
    fill = NA, label.size = 0, label.r = unit(0, "pt"),
    lineheight = 1.05
  ) +

  coord_cartesian(xlim = range(d_multi$ln_genome_size_c, na.rm = TRUE), clip = "off") +
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = c(0.00, 0.25, 0.50, 0.75, 1.00),
    labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +

  theme_classic(base_size = 16) +
  theme(
    text = element_text(family = "Arial", colour = "black"),
    axis.title.x = element_text(size = 18, margin = margin(t = 8)),
    axis.title.y = element_text(size = 18, margin = margin(r = 10)),
    axis.text = element_text(size = 14),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.4),
    axis.ticks.length = unit(5, "pt"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.3),
    plot.margin = margin(10, 18, 10, 28),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.key.width = unit(22, "pt"),
    legend.box.margin = margin(t = 4, r = 0, b = 0, l = 0)
  ) +
  guides(linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(
    x = "ln(Genome size [bp]) (centered)",
    y = "Intergenic proportion",
    linetype = NULL
  )

ggsave("SupFig_intergenic_prop_LS_GS_interaction.png",
       pC, width = 6.2, height = 6.6, dpi = 300, limitsize = FALSE)

cat("Saved: SupFig_intergenic_prop_LS_GS_interaction.png\n")

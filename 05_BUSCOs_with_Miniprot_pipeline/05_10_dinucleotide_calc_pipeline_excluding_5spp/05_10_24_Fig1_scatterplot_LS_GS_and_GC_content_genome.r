## ============================================================
## Main Fig 1 (GENOME ONLY): three scatterplots with PGLS stats
##   Fig1a: ln(Genome size [bp]) vs ln(Lifespan [days])
##   Fig1b: ln(Genome size [bp]) vs Genome GC content
##   Fig1c: ln(Lifespan [days]) vs Genome GC content
##
## - genome rows only (Region == "genome")
## - GC content = G_Rate + C_Rate
## - Each panel shows: PGLS slope (β), p, lambda (λ), n
## - PGLS line shown (solid)
## - Arial, black
## - Output: single PNG with 3 panels horizontally
## ============================================================

library(tidyverse)
library(grid)       # unit()
library(ape)
library(nlme)
library(patchwork)

## ---- INPUT ----
input_file <- "filtered_species_with_tip_name_dinuc_Github_check_plus_intergenic.csv"
df <- read.csv(input_file, check.names = FALSE)

## ---- CHECK FORMAT ----
required_cols <- c("Region", "Tip_Name", "Average_lifespan_days", "Genome_size_bp", "G_Rate", "C_Rate")
missing <- setdiff(required_cols, colnames(df))
if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

## ---- FILTER: GENOME ONLY ----
df <- df %>%
  mutate(Region = as.character(Region)) %>%
  filter(Region == "genome")

if (nrow(df) == 0) stop("No genome rows found (Region == 'genome').")

## ---- DERIVED VARS ----
df <- df %>%
  mutate(
    species_name = as.character(Tip_Name),

    Average_lifespan_days = suppressWarnings(as.numeric(Average_lifespan_days)),
    Genome_size_bp        = suppressWarnings(as.numeric(Genome_size_bp)),
    G_Rate                = suppressWarnings(as.numeric(G_Rate)),
    C_Rate                = suppressWarnings(as.numeric(C_Rate)),

    ln_lifespan = ifelse(!is.na(Average_lifespan_days) & Average_lifespan_days > 0,
                         log(Average_lifespan_days), NA_real_),
    ln_genome_size = ifelse(!is.na(Genome_size_bp) & Genome_size_bp > 0,
                            log(Genome_size_bp), NA_real_),

    GC_content = G_Rate + C_Rate
  )

if (all(is.na(df$GC_content))) stop("GC_content all NA. Check G_Rate and C_Rate.")
if (any(df$GC_content < 0 | df$GC_content > 1, na.rm = TRUE)) {
  warning("Some GC_content values outside [0,1]. Check if rates are normalized.")
}

## ---- LOAD TREE ----
phylo_tree <- read.tree("/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_exon_intron_separated/gene_region_separated_sequences_no_internal_stop_3divisible/tetranucleotide_count/DNA_mechanics/PGLS/labelled_supertree_ottnames.tre")
phylo_tree <- compute.brlen(phylo_tree, method = "Grafen")

## ---- MATCH TO TREE ----
df <- df %>% filter(!is.na(species_name), species_name %in% phylo_tree$tip.label)
if (nrow(df) < 3) stop("Too few genome species after matching to tree (n < 3).")

phylo_tree <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, df$species_name))

## ============================================================
## Helper: Fit PGLS y ~ x, return slope, p, lambda, n, intercept
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
    return(list(
      df = d, n = n,
      intercept = NA_real_, slope = NA_real_, p = NA_real_, lambda = NA_real_
    ))
  }

  tr <- drop.tip(tree, setdiff(tree$tip.label, d$species_name))

  fit <- gls(
    y ~ x,
    data = d,
    correlation = corPagel(1, phy = tr, form = ~species_name),
    method = "ML"
  )

  tt <- summary(fit)$tTable
  slope <- tt["x", "Value"]
  pval  <- tt["x", "p-value"]
  intercept <- tt["(Intercept)", "Value"]

  lambda_value <- tryCatch(
    as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE)),
    error = function(e) NA_real_
  )

  list(df = d, n = n, intercept = intercept, slope = slope, p = pval, lambda = lambda_value)
}

## ============================================================
## Helper: plot a panel with PGLS line + annotation
## ============================================================
plot_panel <- function(fit_obj, xlab, ylab,
                       y_is_gc = FALSE,
                       y_fixed = NULL) {


  d <- fit_obj$df

  ann <- paste0(
    "PGLS slope = ", formatC(fit_obj$slope, format = "f", digits = 3), "\n",
    " p = ", formatC(fit_obj$p, format = "e", digits = 2), "\n",
    "  λ = ", formatC(fit_obj$lambda, format = "f", digits = 2), "\n",
    "  n = ", fit_obj$n
  )

  # X range: add 5% padding
  x_rng <- range(d$x, na.rm = TRUE)
  x_pad <- diff(x_rng) * 0.05
  x_limits <- c(x_rng[1] - x_pad, x_rng[2] + x_pad)

  # Y range
  if (!is.null(y_fixed)) {
    y_limits <- y_fixed
  } else if (y_is_gc) {
    y_limits <- c(0, 1)
  } else {
    y_rng <- range(d$y, na.rm = TRUE)
    y_pad <- diff(y_rng) * 0.05
    y_limits <- c(y_rng[1] - y_pad, y_rng[2] + y_pad)
  }


  ggplot(d, aes(x = x, y = y)) +
    geom_point(size = 1.2, alpha = 0.50, color = "black") +

    # PGLS regression line
    geom_abline(
      intercept = fit_obj$intercept,
      slope = fit_obj$slope,
      colour = "black",
      linewidth = 1.2, alpha = 0.5
    ) +

    # Annotation (top-left)
    geom_text(
      data = tibble(x = -Inf, y = Inf, label = ann),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = -0.1,
      vjust = 1.2,
      size = 4.0,
      color = "black",
      lineheight = 1.05
    ) +

    coord_cartesian(xlim = x_limits, ylim = y_limits) +

    theme_classic(base_size = 16) +
    theme(
      text = element_text(family = "Arial", colour = "black"),
      axis.title.x = element_text(size = 18, margin = margin(t = 8)),
      axis.title.y = element_text(size = 18, margin = margin(r = 8)),
      axis.text = element_text(size = 14),

      # Remove axis lines
      axis.line = element_blank(),

      axis.ticks = element_line(linewidth = 0.4),
      axis.ticks.length = unit(5, "pt"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.3),
      plot.margin = margin(10, 18, 10, 18)
    ) +
    labs(x = xlab, y = ylab)
}


## ============================================================
## Fit models + make panels
## ============================================================

# Fig1a: ln(genome size) vs ln(lifespan)
fit_a <- fit_pgls(df, x = "ln_genome_size", y = "ln_lifespan", tree = phylo_tree)
p1a <- plot_panel(
  fit_a,
  xlab = "ln(Genome size [bp])",
  ylab = "ln(Lifespan [days])",
  y_fixed = c(0, 15)
)

# Fig1b: ln(genome size) vs genome GC content
fit_b <- fit_pgls(df, x = "ln_genome_size", y = "GC_content", tree = phylo_tree)
p1b <- plot_panel(
  fit_b,
  xlab = "ln(Genome size [bp])",
  ylab = "GC content",
  y_is_gc = TRUE
)

# Fig1c: ln(lifespan) vs genome GC content
fit_c <- fit_pgls(df, x = "ln_lifespan", y = "GC_content", tree = phylo_tree)
p1c <- plot_panel(
  fit_c,
  xlab = "ln(Lifespan [days])",
  ylab = "GC content",
  y_is_gc = TRUE
)

## ============================================================
## Combine + save
## ============================================================
p_combined <- p1a | p1b | p1c

ggsave(
  "Fig1_genome_lnGenome_lnLifespan_GCcontent_PGLS.png",
  p_combined,
  width = 18.5,
  height = 6.0,
  dpi = 300
)

cat("Saved: Fig1_genome_lnGenome_lnLifespan_GCcontent_PGLS.png\n")

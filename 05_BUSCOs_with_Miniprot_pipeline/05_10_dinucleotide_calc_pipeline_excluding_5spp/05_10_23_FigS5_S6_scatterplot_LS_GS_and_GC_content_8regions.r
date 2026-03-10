## ============================================================
## Scatterplots (8 regions): GC content vs
##   (A) ln(lifespan [days])
##   (B) ln(genome size [bp])
##
## - Per-panel PGLS (univariate): GC_content ~ x
## - Show PGLS regression line (solid) + annotate β and FDR p (moved inward)
## - y label: "GC content"
##
## Input: filtered_species_with_tip_name_dinuc_Github_check_plus_intergenic.csv
## GC content = G_Rate + C_Rate
## Output: 2 PNGs
## ============================================================

library(tidyverse)
library(grid)   # unit()
library(ape)
library(nlme)

## ---- INPUT ----
input_file <- "filtered_species_with_tip_name_dinuc_Github_check_plus_intergenic.csv"
df <- read.csv(input_file, check.names = FALSE)

## ---- CHECK FORMAT ----
required_cols <- c("Region", "Average_lifespan_days", "Genome_size_bp", "G_Rate", "C_Rate", "Tip_Name")
missing <- setdiff(required_cols, colnames(df))
if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

## ---- REGIONS (8) ----
region_levels <- c("genome", "intergenic", "upstream1", "upstream2",
                   "exon", "intron", "downstream1", "downstream2")
df$Region <- as.character(df$Region)
df <- df %>% filter(Region %in% region_levels)

missing_regions <- setdiff(region_levels, unique(df$Region))
if (length(missing_regions) > 0) {
  stop("These requested regions are missing from the input file: ",
       paste(missing_regions, collapse = ", "),
       "\nPresent regions are: ",
       paste(sort(unique(df$Region)), collapse = ", "))
}

## ---- NUMERIC + DERIVED VARS ----
df <- df %>%
  mutate(
    Average_lifespan_days = suppressWarnings(as.numeric(Average_lifespan_days)),
    Genome_size_bp        = suppressWarnings(as.numeric(Genome_size_bp)),
    G_Rate                = suppressWarnings(as.numeric(G_Rate)),
    C_Rate                = suppressWarnings(as.numeric(C_Rate)),

    ln_lifespan = ifelse(!is.na(Average_lifespan_days) & Average_lifespan_days > 0,
                         log(Average_lifespan_days), NA_real_),
    ln_genome_size = ifelse(!is.na(Genome_size_bp) & Genome_size_bp > 0,
                            log(Genome_size_bp), NA_real_),

    GC_content = G_Rate + C_Rate,

    species_name = as.character(Tip_Name),
    Region = factor(Region, levels = region_levels)
  )

if (all(is.na(df$GC_content))) stop("GC_content is all NA. Check G_Rate and C_Rate values.")
if (any(df$GC_content < 0 | df$GC_content > 1, na.rm = TRUE)) {
  warning("Some GC_content values are outside [0,1]. Check if rates are correctly normalized.")
}

## ---- LOAD TREE ----
phylo_tree <- read.tree("/home/y-kondo/protostome_lifespan_model/BUSCO_seq_directory/up_downstreams_exon_intron_separated/gene_region_separated_sequences_no_internal_stop_3divisible/tetranucleotide_count/DNA_mechanics/PGLS/labelled_supertree_ottnames.tre")
phylo_tree <- compute.brlen(phylo_tree, method = "Grafen")

## ---- LABEL HELPERS ----
facet_title <- function(x) {
  x <- as.character(x)
  x <- gsub("([a-zA-Z]+)([0-9]+)", "\\1 \\2", x)  # upstream1 -> upstream 1
  tools::toTitleCase(x)
}

## ============================================================
## PGLS helper: per-region univariate model GC_content ~ x
## Returns Region, beta, p, n, FDR p, label
## ============================================================
run_pgls_by_region <- function(dat, x_col, phylo_tree) {
  res <- list()

  for (r in levels(dat$Region)) {
    d_r <- dat %>%
      filter(Region == r) %>%
      transmute(
        species_name = as.character(species_name),
        GC_content   = as.numeric(GC_content),
        x            = as.numeric(.data[[x_col]])
      ) %>%
      drop_na(species_name, GC_content, x) %>%
      filter(is.finite(GC_content), is.finite(x)) %>%
      filter(species_name %in% phylo_tree$tip.label)

    if (nrow(d_r) < 3) {
      res[[r]] <- tibble(Region = r, beta = NA_real_, p = NA_real_, n = nrow(d_r))
      next
    }

    tr <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, d_r$species_name))

    fit <- tryCatch(
      gls(
        GC_content ~ x,
        data = d_r,
        correlation = corPagel(1, phy = tr, form = ~species_name),
        method = "ML"
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      res[[r]] <- tibble(Region = r, beta = NA_real_, p = NA_real_, n = nrow(d_r))
      next
    }

    tt <- summary(fit)$tTable
    beta <- if ("x" %in% rownames(tt)) tt["x", "Value"] else NA_real_
    pval <- if ("x" %in% rownames(tt)) tt["x", "p-value"] else NA_real_

    res[[r]] <- tibble(Region = r, beta = beta, p = pval, n = nrow(d_r))
  }

  out <- bind_rows(res) %>%
    mutate(
      Region = factor(Region, levels = levels(dat$Region)),
      p_fdr = p.adjust(p, method = "BH"),
      label = paste0(
        "β = ", formatC(beta, format = "f", digits = 3), "\n",
        "FDR p = ", formatC(p_fdr, format = "e", digits = 2)
      )
    )

  out
}

## ============================================================
## PGLS line helper: get intercept + slope per region
## ============================================================
fit_pgls_lines_by_region <- function(dat, x_col, phylo_tree) {
  res <- list()

  for (r in levels(dat$Region)) {
    d_r <- dat %>%
      filter(Region == r) %>%
      transmute(
        species_name = as.character(species_name),
        GC_content   = as.numeric(GC_content),
        x            = as.numeric(.data[[x_col]])
      ) %>%
      drop_na(species_name, GC_content, x) %>%
      filter(is.finite(GC_content), is.finite(x)) %>%
      filter(species_name %in% phylo_tree$tip.label)

    if (nrow(d_r) < 3) next

    tr <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, d_r$species_name))

    fit <- tryCatch(
      gls(
        GC_content ~ x,
        data = d_r,
        correlation = corPagel(1, phy = tr, form = ~species_name),
        method = "ML"
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    cf <- coef(fit)
    res[[r]] <- tibble(
      Region = r,
      intercept = unname(cf["(Intercept)"]),
      slope = unname(cf["x"])
    )
  }

  bind_rows(res) %>%
    mutate(Region = factor(Region, levels = levels(dat$Region)))
}

## ============================================================
## Plot function with PGLS solid line + inward annotation
## ============================================================
make_gc_plot_with_pgls <- function(dat, x_col, x_lab, out_png, phylo_tree) {

  d <- dat %>%
    transmute(
      Region = Region,
      species_name = species_name,
      GC_content = as.numeric(GC_content),
      x = as.numeric(.data[[x_col]])
    ) %>%
    drop_na(species_name, x, GC_content) %>%
    filter(is.finite(x), is.finite(GC_content))

  if (nrow(d) < 10) stop("Too few complete rows to plot for: ", x_col)

  # PGLS stats (β + FDR p) per region
  stats_df <- run_pgls_by_region(dat, x_col, phylo_tree)

  # PGLS line parameters per region
  line_df <- fit_pgls_lines_by_region(dat, x_col, phylo_tree)

  p <- ggplot(d, aes(x = x, y = GC_content)) +
    geom_point(size = 0.9, alpha = 0.45, color = "black") +

    # PGLS regression line (solid)
    geom_abline(
      data = line_df,
      aes(intercept = intercept, slope = slope),
      colour = "black",
      linewidth = 1.2
    ) +

    facet_wrap(~Region, ncol = 4, labeller = labeller(Region = facet_title)) +

    # PGLS annotation moved inward
    geom_text(
      data = stats_df,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.1,
      vjust = 1.2,
      size = 4.0,
      color = "black",
      lineheight = 1.05
    ) +

    scale_x_continuous(expand = expansion(mult = c(0.06, 0.06))) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.06)), limits = c(0, 1)) +

    theme_classic(base_size = 18) +
    theme(
      text = element_text(family = "Arial", colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(size = 16, face = "bold", colour = "black"),
      axis.title.x = element_text(size = 22, margin = margin(t = 10), colour = "black"),
      axis.title.y = element_text(size = 22, margin = margin(r = 10), colour = "black"),
      axis.text = element_text(size = 14, colour = "black"),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks.length = unit(5, "pt"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.0),
      plot.margin = margin(8, 8, 8, 8)
    ) +
    labs(
      x = x_lab,
      y = "GC content"
    )

  ggsave(out_png, p, width = 16, height = 8.5, dpi = 300)
  cat("Saved: ", out_png, "\n", sep = "")
}

## ---- (A) ln(lifespan) vs GC ----
make_gc_plot_with_pgls(
  dat = df,
  x_col = "ln_lifespan",
  x_lab = "ln(Lifespan [days])",
  out_png = "Scatter_GCcontent_vs_lnLifespan_8regions_PGLSline.png",
  phylo_tree = phylo_tree
)

## ---- (B) ln(genome size) vs GC ----
make_gc_plot_with_pgls(
  dat = df,
  x_col = "ln_genome_size",
  x_lab = "ln(Genome size [bp])",
  out_png = "Scatter_GCcontent_vs_lnGenomeSize_8regions_PGLSline.png",
  phylo_tree = phylo_tree
)

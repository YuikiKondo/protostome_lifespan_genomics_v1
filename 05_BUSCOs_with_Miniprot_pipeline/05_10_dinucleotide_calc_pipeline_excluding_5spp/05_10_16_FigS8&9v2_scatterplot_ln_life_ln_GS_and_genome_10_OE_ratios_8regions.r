## ============================================================
## Supplementary scatterplots (ALL 8 regions), batched
##   (A) ln(lifespan) vs dinuc O/E
##   (B) ln(genome size [bp]) vs dinuc O/E
##
## Adds: OLS equation + R^2 in top-left of each panel (NO p)
## ============================================================

library(tidyverse)
library(grid)  # unit()

## ---- INPUT ----
input_file <- "di_obs_exp_ratio_rev_comp_collapsed_only_with_lifespan_tip_names_and_genome_size_Github_check.csv"
data <- read.csv(input_file, check.names = FALSE)

## ---- BASIC PREP ----
data$Average_lifespan_days <- suppressWarnings(as.numeric(data$Average_lifespan_days))
data$ln_lifespan <- ifelse(!is.na(data$Average_lifespan_days) & data$Average_lifespan_days > 0,
                           log(data$Average_lifespan_days), NA_real_)

data$Genome_size_bp <- suppressWarnings(as.numeric(data$Genome_size_bp))
data$ln_genome_size <- ifelse(!is.na(data$Genome_size_bp) & data$Genome_size_bp > 0,
                              log(data$Genome_size_bp), NA_real_)

data$Region <- as.character(data$Region)

## ---- TARGET REGIONS ----
region_levels <- c("genome", "intergenic", "upstream1", "upstream2", "exon", "intron", "downstream1", "downstream2")
data <- data %>% filter(Region %in% region_levels)

missing_regions <- setdiff(region_levels, unique(data$Region))
if (length(missing_regions) > 0) {
  stop("These requested regions are missing from the input file: ",
       paste(missing_regions, collapse = ", "))
}

## ---- MOTIF LEVELS (16) ----
motif_levels <- c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")

## ---- COLLAPSED → EXPAND MAP (ONLY for genome + intergenic) ----
collapse_map <- list(
  "AT_OE" = c("AT"),
  "TA_OE" = c("TA"),
  "CG_OE" = c("CG"),
  "GC_OE" = c("GC"),
  "CA_and_TG_OE" = c("CA","TG"),
  "AC_and_GT_OE" = c("AC","GT"),
  "AG_and_CT_OE" = c("AG","CT"),
  "GA_and_TC_OE" = c("GA","TC"),
  "AA_and_TT_OE" = c("AA","TT"),
  "CC_and_GG_OE" = c("CC","GG")
)

collapsed_regions <- c("genome", "intergenic")

## ---- Expected OE column sets ----
bases <- c("A","C","G","T")
oe16_cols <- as.vector(outer(bases, bases, paste0))
oe16_cols <- paste0(oe16_cols, "_OE")  # AA_OE ... TT_OE
collapsed10_cols <- names(collapse_map)

## ---- Build long table: species × region × motif ----
## Collapsed regions: expand 10 -> 16 by duplication
collapsed_long <- data %>%
  filter(Region %in% collapsed_regions) %>%
  select(Region, ln_lifespan, ln_genome_size, all_of(intersect(collapsed10_cols, colnames(data)))) %>%
  pivot_longer(
    cols = all_of(intersect(collapsed10_cols, colnames(data))),
    names_to = "CollapsedCol",
    values_to = "OE"
  ) %>%
  mutate(OE = suppressWarnings(as.numeric(OE))) %>%
  drop_na(OE) %>%
  rowwise() %>%
  do({
    members <- collapse_map[[.$CollapsedCol]]
    tibble(
      Region = .$Region,
      ln_lifespan = .$ln_lifespan,
      ln_genome_size = .$ln_genome_size,
      Motif = members,
      OE = .$OE
    )
  }) %>%
  ungroup()

## Non-collapsed regions: use 16 motif-specific columns
noncollapsed_long <- data %>%
  filter(!Region %in% collapsed_regions) %>%
  select(Region, ln_lifespan, ln_genome_size, all_of(intersect(oe16_cols, colnames(data)))) %>%
  pivot_longer(
    cols = all_of(intersect(oe16_cols, colnames(data))),
    names_to = "MotifCol",
    values_to = "OE"
  ) %>%
  mutate(
    OE = suppressWarnings(as.numeric(OE)),
    Motif = gsub("_OE$", "", MotifCol)
  ) %>%
  select(Region, ln_lifespan, ln_genome_size, Motif, OE) %>%
  drop_na(OE)

plot_long <- bind_rows(collapsed_long, noncollapsed_long) %>%
  mutate(
    Region = factor(Region, levels = region_levels),
    Motif  = factor(Motif, levels = motif_levels)
  ) %>%
  filter(!is.na(Region), !is.na(Motif))

## ---- Motif batches: 8,8 ----
motif_batches <- list(
  batch1 = motif_levels[1:8],
  batch2 = motif_levels[9:16]
)

## ---- Output dirs ----
dir.create("SuppFig_batched_lnLifespan_vs_OE", showWarnings = FALSE)
dir.create("SuppFig_batched_lnGenomeSize_vs_OE", showWarnings = FALSE)

## ------------------------------------------------------------
## NEW: build labels per facet (Motif × Region)
## ------------------------------------------------------------
make_lm_labels_grid <- function(df, motifs, x_col, y_col = "OE", min_n = 8) {

  df_use <- df %>%
    filter(Motif %in% motifs) %>%
    filter(!is.na(.data[[x_col]]), is.finite(.data[[x_col]])) %>%
    filter(!is.na(.data[[y_col]]), is.finite(.data[[y_col]]))

  lab_df <- df_use %>%
    group_by(Motif, Region) %>%
    summarise(
      n = n(),
      .groups = "drop_last"
    ) %>%
    mutate(ok = n >= min_n) %>%
    select(-n) %>%
    right_join(df_use %>% distinct(Motif, Region), by = c("Motif", "Region")) %>%
    mutate(ok = ifelse(is.na(ok), FALSE, ok)) %>%
    group_by(Motif, Region) %>%
    group_modify(~{
      if (!.x$ok[1]) {
        return(tibble(label = ""))
      }
      fit <- lm(reformulate(x_col, y_col), data = df_use %>% filter(Motif == .y$Motif, Region == .y$Region))
      sm  <- summary(fit)

      a  <- unname(coef(fit)[1])
      b  <- unname(coef(fit)[2])
      r2 <- unname(sm$r.squared)

      label <- sprintf("y = %.3f %s %.3f x\nR² = %.3f",
                       a,
                       ifelse(b >= 0, "+", "−"),
                       abs(b),
                       r2)
      tibble(label = label)
    }) %>%
    ungroup() %>%
    mutate(
      Motif  = factor(Motif, levels = motifs),
      Region = factor(Region, levels = region_levels)
    )

  lab_df
}

## ------------------------------------------------------------
## Plot helper: one batch -> motif rows × region columns
## - Region strip labels: TitleCase
## - Motif strip labels: "TT O/E"
## - Add equation + R^2 in top-left of each panel
## ------------------------------------------------------------
make_batch_plot <- function(df, motifs, x_col, x_lab, out_png) {

  d <- df %>%
    filter(Motif %in% motifs) %>%
    filter(!is.na(.data[[x_col]]), is.finite(.data[[x_col]])) %>%
    drop_na(OE) %>%
    filter(is.finite(OE)) %>%
    mutate(
      Motif = factor(Motif, levels = motifs),   # keep batch order
      Region = factor(Region, levels = region_levels)
    )

  if (nrow(d) < 20) {
    message("Skipping ", out_png, " (too few points overall: n=", nrow(d), ")")
    return(invisible(NULL))
  }

  # Use a single y-limit across THIS BATCH
  y_lim <- range(d$OE, na.rm = TRUE)

  # Labels per facet (Motif × Region)
  ann_df <- make_lm_labels_grid(d, motifs = motifs, x_col = x_col, y_col = "OE", min_n = 8)

  motif_labeller  <- function(x) paste0(x, " O/E")          # "TT O/E"
  region_labeller <- function(x) tools::toTitleCase(x)      # "Genome", "Upstream1", ...

  p <- ggplot(d, aes(x = .data[[x_col]], y = OE)) +
    geom_point(size = 0.7, alpha = 0.40) +
    geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.8, alpha = 0.8) +
    geom_text(
      data = ann_df,
      aes(label = label),
      x = -Inf,
      y = Inf,
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1.15,
      size = 4.2,
      family = "Arial",
      colour = "black",
      lineheight = 0.95
    ) +
    facet_grid(
      rows = vars(Motif),
      cols = vars(Region),
      labeller = labeller(
        Motif = motif_labeller,
        Region = region_labeller
      )
    ) +
    coord_cartesian(ylim = y_lim, clip = "off") +
    scale_x_continuous(expand = expansion(mult = c(0.07, 0.07))) +
    scale_y_continuous(expand = expansion(mult = c(0.07, 0.07))) +
    theme_classic(base_size = 18) +
    theme(
      text = element_text(family = "Arial", colour = "black"),

      strip.background = element_blank(),
      strip.text.x = element_text(size = 20, face = "bold", colour = "black"),
      strip.text.y = element_text(size = 20, face = "bold", colour = "black"),

      axis.title.x = element_text(size = 28, margin = margin(t = 10), colour = "black"),
      axis.title.y = element_text(size = 28, margin = margin(r = 10), colour = "black"),
      axis.text = element_text(size = 13, colour = "black"),

      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.length = unit(4, "pt"),

      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
      plot.margin = margin(8, 8, 8, 8)
    ) +
    labs(
      x = x_lab,
      y = "Dinucleotide O/E ratio"
    )

  # Height scales with number of motif rows
  n_rows <- length(motifs)
  out_h <- 2.5 + 2.5 * n_rows
  out_w <- 22

  ggsave(out_png, p, width = out_w, height = out_h, dpi = 300)
  cat("Saved: ", out_png, "\n", sep = "")
}

## ------------------------------------------------------------
## Generate batched figures
## ------------------------------------------------------------
for (bn in names(motif_batches)) {
  motifs <- motif_batches[[bn]]

  ## (A) ln(lifespan)
  make_batch_plot(
    df = plot_long,
    motifs = motifs,
    x_col = "ln_lifespan",
    x_lab = "ln(Lifespan [days])",
    out_png = file.path("SuppFig_batched_lnLifespan_vs_OE",
                        paste0("SuppFig_lnLifespan_vs_OE_", bn, "_", length(motifs), "motifs.png"))
  )

  ## (B) ln(genome size)
  make_batch_plot(
    df = plot_long,
    motifs = motifs,
    x_col = "ln_genome_size",
    x_lab = "ln(Genome size [bp])",
    out_png = file.path("SuppFig_batched_lnGenomeSize_vs_OE",
                        paste0("SuppFig_lnGenomeSize_vs_OE_", bn, "_", length(motifs), "motifs.png"))
  )
}

cat("Done.\nOutputs:\n",
    " - SuppFig_batched_lnLifespan_vs_OE/ (2 PNGs: 8,8 motifs)\n",
    " - SuppFig_batched_lnGenomeSize_vs_OE/ (2 PNGs: 8,8 motifs)\n",
    sep = "")
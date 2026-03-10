## ============================================================
## Supplementary scatterplots (ALL 8 regions), batched
##   raw dinucleotide frequency (Rate) vs dinuc O/E
##
## IMPORTANT:
## - This input already has 16 motif columns for ALL regions.
## - For genome/intergenic, reverse-complement collapsed values are duplicated
##   (AA == TT, AC == GT, etc.), so those panels will be identical automatically.
##
## Add: Pearson's r and BH-FDR p-value per panel
## FDR is corrected globally across all 128 panels.
## ============================================================

library(tidyverse)
library(grid)  # unit()

## ---- INPUT ----
input_file <- "filtered_species_with_tip_name_dinuc_Github_check_plus_intergenic.csv"
data <- read.csv(input_file, check.names = FALSE)

## ---- REMOVE MONONUCLEOTIDE RATE COLUMNS (just in case) ----
mono_rate_cols <- c("A_Rate", "T_Rate", "G_Rate", "C_Rate")
data <- data %>% select(-any_of(mono_rate_cols))

data$Region <- as.character(data$Region)

## ---- TARGET REGIONS ----
region_levels <- c(
  "genome", "intergenic",
  "upstream1", "upstream2",
  "exon", "intron",
  "downstream1", "downstream2"
)
data <- data %>% filter(Region %in% region_levels)

missing_regions <- setdiff(region_levels, unique(data$Region))
if (length(missing_regions) > 0) {
  stop("These requested regions are missing from the input file: ",
       paste(missing_regions, collapse = ", "))
}

## ---- MOTIF LEVELS (16) ----
motif_levels <- c(
  "AA","AC","AG","AT","CA","CC","CG","CT",
  "GA","GC","GG","GT","TA","TC","TG","TT"
)

## ---- Expected column sets (16 motifs) ----
bases <- c("A","C","G","T")
rate16_cols <- paste0(as.vector(outer(bases, bases, paste0)), "_Rate")
oe16_cols   <- paste0(as.vector(outer(bases, bases, paste0)), "_OE")

## Sanity check: ensure columns exist
missing_rate <- setdiff(rate16_cols, colnames(data))
missing_oe   <- setdiff(oe16_cols,   colnames(data))
if (length(missing_rate) > 0 || length(missing_oe) > 0) {
  stop("Missing motif columns.\n",
       "Missing *_Rate: ", paste(missing_rate, collapse = ", "), "\n",
       "Missing *_OE: ", paste(missing_oe, collapse = ", "))
}

## Choose a per-species ID column (not strictly needed for plotting)
id_col <- if ("Tip_Name" %in% colnames(data)) "Tip_Name" else "Organism Name"
if (!id_col %in% colnames(data)) stop("Neither Tip_Name nor Organism Name exists in the file.")

## ============================================================
## Build long table: region × motif with Rate and OE
## (No special handling for genome/intergenic; duplication already present)
## ============================================================
plot_long <- data %>%
  select(Region, all_of(id_col), all_of(rate16_cols), all_of(oe16_cols)) %>%
  pivot_longer(
    cols = -c(Region, all_of(id_col)),
    names_to = "Col",
    values_to = "Value"
  ) %>%
  mutate(
    Value  = suppressWarnings(as.numeric(Value)),
    Suffix = ifelse(grepl("_Rate$", Col), "Rate", "OE"),
    Motif  = gsub("_(Rate|OE)$", "", Col)
  ) %>%
  select(Region, all_of(id_col), Motif, Suffix, Value) %>%
  pivot_wider(
    id_cols = c(Region, all_of(id_col), Motif),
    names_from = Suffix,
    values_from = Value,
    values_fn = list(Value = mean) # defensive
  ) %>%
  drop_na(Rate, OE) %>%
  mutate(
    Rate = as.numeric(Rate),
    OE   = as.numeric(OE),
    Region = factor(Region, levels = region_levels),
    Motif  = factor(Motif,  levels = motif_levels)
  ) %>%
  filter(is.finite(Rate), is.finite(OE), !is.na(Region), !is.na(Motif))

## ---- Motif batches: 8,8 ----
motif_batches <- list(
  batch1 = motif_levels[1:8],
  batch2 = motif_levels[9:16]
)

## ---- Output dir ----
dir.create("SuppFig_batched_Rate_vs_OE", showWarnings = FALSE)

## ------------------------------------------------------------
## Compute Pearson r + p for ALL panels, then BH-FDR globally
## ------------------------------------------------------------
panel_stats_all <- plot_long %>%
  group_by(Motif, Region) %>%
  summarise(
    n = sum(is.finite(Rate) & is.finite(OE)),
    r = suppressWarnings(cor(Rate, OE, method = "pearson")),
    p = suppressWarnings(cor.test(Rate, OE, method = "pearson")$p.value),
    .groups = "drop"
  ) %>%
  mutate(
    p_fdr = p.adjust(p, method = "BH"),
    label = sprintf("r = %.3f\nFDR p = %.3g", r, p_fdr)
  )

## ------------------------------------------------------------
## Plot helper: one batch
## ------------------------------------------------------------
make_batch_plot <- function(df, motifs, out_png) {

  d <- df %>%
    filter(Motif %in% motifs) %>%
    mutate(
      Motif  = factor(Motif, levels = motifs),
      Region = factor(Region, levels = region_levels)
    ) %>%
    drop_na(Rate, OE) %>%
    filter(is.finite(Rate), is.finite(OE))

  if (nrow(d) < 20) {
    message("Skipping ", out_png, " (too few points)")
    return(invisible(NULL))
  }

  # per-panel labels subset for this batch
  ann_df <- panel_stats_all %>%
    filter(Motif %in% motifs) %>%
    mutate(
      Motif  = factor(Motif, levels = motifs),
      Region = factor(Region, levels = region_levels)
    )

  y_lim <- range(d$OE, na.rm = TRUE)

  p <- ggplot(d, aes(x = Rate, y = OE)) +
    geom_point(size = 0.7, alpha = 0.40) +
    geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.8) +
    geom_text(
      data = ann_df,
      aes(label = label),
      x = -Inf, y = Inf,
      inherit.aes = FALSE,
      hjust = -0.05, vjust = 1.15,
      size = 3.6,
      family = "Arial",
      colour = "black",
      lineheight = 0.95
    ) +
    facet_grid(
      rows = vars(Motif),
      cols = vars(Region),
      labeller = labeller(
        Motif  = function(x) paste0(x, " O/E"),
        Region = function(x) tools::toTitleCase(x)
      )
    ) +
    coord_cartesian(ylim = y_lim, clip = "off") +
    scale_x_continuous(expand = expansion(mult = c(0.07, 0.07))) +
    scale_y_continuous(expand = expansion(mult = c(0.07, 0.07))) +
    theme_classic(base_size = 18) +
    theme(
      text = element_text(family = "Arial", colour = "black"),
      strip.background = element_blank(),
      strip.text.x = element_text(size = 20, face = "bold"),
      strip.text.y = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 28),
      axis.title.y = element_text(size = 28),
      axis.text = element_text(size = 13),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
      plot.margin = margin(8, 8, 8, 8)
    ) +
    labs(
      x = "Raw dinucleotide frequency (Rate)",
      y = "Dinucleotide O/E ratio"
    )

  out_h <- 2.5 + 2.5 * length(motifs)
  ggsave(out_png, p, width = 22, height = out_h, dpi = 300)
  cat("Saved:", out_png, "\n")
}

## ------------------------------------------------------------
## Generate figures
## ------------------------------------------------------------
for (bn in names(motif_batches)) {
  make_batch_plot(
    plot_long,
    motif_batches[[bn]],
    file.path("SuppFig_batched_Rate_vs_OE",
              paste0("SuppFig_Rate_vs_OE_", bn, "_", length(motif_batches[[bn]]), "motifs.png"))
  )
}

cat("Done.\nOutputs: SuppFig_batched_Rate_vs_OE/\n")


cat(
  "Significant (FDR < 0.05):",
  sum(panel_stats_all$p_fdr < 0.05, na.rm = TRUE),
  "/",
  nrow(panel_stats_all),
  "\n"
)


panel_stats_all %>%
  filter(p_fdr > 0.05) %>%
  arrange(Region, Motif)
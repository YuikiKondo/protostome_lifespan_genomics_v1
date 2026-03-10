## ============================================================
## 3 heatmaps (16 motifs × 6 regions) from PGLS (ln(O/E) response)
## - genome + intergenic: collapsed 10 motifs expanded to 16
## - Regions: Genome, Intergenic, Upstream, Exon, Intron, Downstream
## - Panels:
##   β1 lifespan, β2 genome size, β3 interaction
## - Each heatmap has its OWN legend on the right
## - All text: Arial, black
## ============================================================

library(tidyverse)
library(grid)
library(patchwork)

## INPUT ----
in_file <- "pgls_results_8_regions_genome10_intergenic10_others16_dinuc_ln_OE_LS_CENTERED_G_CENTERED_interaction_OE_as_response.csv"
df <- read.csv(in_file, check.names = FALSE)

## CHECK REQUIRED COLS ----
required <- c(
  "Region", "Dinucleotide",
  "Coefficient_ln_lifespan_c", "p.adjusted.fdr_ln_lifespan_c",
  "Coefficient_ln_Genome_Size_c", "p.adjusted.fdr_ln_Genome_Size_c",
  "Coefficient_Interaction", "p.adjusted.fdr_Interaction"
)
missing <- setdiff(required, colnames(df))
if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

## CLEAN + COERCE ----
df <- df %>%
  mutate(
    Region = as.character(Region),
    Dinucleotide = as.character(Dinucleotide),
    Dinuc2 = trimws(gsub("_OE$", "", Dinucleotide)),

    Coefficient_ln_lifespan_c = as.numeric(Coefficient_ln_lifespan_c),
    p.adjusted.fdr_ln_lifespan_c = as.numeric(p.adjusted.fdr_ln_lifespan_c),

    Coefficient_ln_Genome_Size_c = as.numeric(Coefficient_ln_Genome_Size_c),
    p.adjusted.fdr_ln_Genome_Size_c = as.numeric(p.adjusted.fdr_ln_Genome_Size_c),

    Coefficient_Interaction = as.numeric(Coefficient_Interaction),
    p.adjusted.fdr_Interaction = as.numeric(p.adjusted.fdr_Interaction)
  )


## REGIONS USED ----
region_levels <- c("genome", "intergenic", "upstream1", "upstream2", "exon", "intron", "downstream1", "downstream2")
df <- df %>% filter(Region %in% region_levels)

## MOTIF LEVELS ----
motif_levels <- c(
  "AA","AC","AG","AT","CA","CC","CG","CT",
  "GA","GC","GG","GT","TA","TC","TG","TT"
)

## ------------------------------------------------------------
## Expand collapsed genome + intergenic
## ------------------------------------------------------------
collapse_map <- list(
  "AT" = "AT", "TA" = "TA", "GC" = "GC", "CG" = "CG",
  "CA_and_TG" = c("CA","TG"),
  "AC_and_GT" = c("AC","GT"),
  "AG_and_CT" = c("AG","CT"),
  "GA_and_TC" = c("GA","TC"),
  "AA_and_TT" = c("AA","TT"),
  "CC_and_GG" = c("CC","GG")
)

collapsed_regions <- c("genome", "intergenic")

expand_rows <- function(df_sub) {
  df_sub %>%
    rowwise() %>%
    do({
      members <- collapse_map[[.$Dinuc2]]
      if (is.null(members)) members <- .$Dinuc2
      tibble(
        Region = .$Region,
        Dinuc2 = members,
        Coefficient_ln_lifespan_c = .$Coefficient_ln_lifespan_c,
        p.adjusted.fdr_ln_lifespan_c = .$p.adjusted.fdr_ln_lifespan_c,
        Coefficient_ln_Genome_Size_c = .$Coefficient_ln_Genome_Size_c,
        p.adjusted.fdr_ln_Genome_Size_c = .$p.adjusted.fdr_ln_Genome_Size_c,
        Coefficient_Interaction = .$Coefficient_Interaction,
        p.adjusted.fdr_Interaction = .$p.adjusted.fdr_Interaction
      )
    }) %>%
    ungroup()
}

collapsed <- df %>% filter(Region %in% collapsed_regions)
non_collapsed <- df %>% filter(!Region %in% collapsed_regions)

all_long <- bind_rows(
  expand_rows(collapsed),
  non_collapsed
) %>%
  mutate(
    Region = factor(Region, levels = region_levels),
    Dinuc2 = factor(Dinuc2, levels = motif_levels)
  ) %>%
  group_by(Region, Dinuc2) %>%
  summarise(
    Coefficient_ln_lifespan_c = mean(Coefficient_ln_lifespan_c, na.rm = TRUE),
    Coefficient_ln_Genome_Size_c = mean(Coefficient_ln_Genome_Size_c, na.rm = TRUE),
    Coefficient_Interaction = mean(Coefficient_Interaction, na.rm = TRUE),
    p.adjusted.fdr_ln_lifespan_c = min(p.adjusted.fdr_ln_lifespan_c, na.rm = TRUE),
    p.adjusted.fdr_ln_Genome_Size_c = min(p.adjusted.fdr_ln_Genome_Size_c, na.rm = TRUE),
    p.adjusted.fdr_Interaction = min(p.adjusted.fdr_Interaction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(Region, Dinuc2)

## ------------------------------------------------------------
## Heatmap function (Arial + black text)
## ------------------------------------------------------------
plot_heatmap <- function(dat, fill_col, sig_col, title_text, legend_title) {

  dat <- dat %>%
    mutate(
      fill_val = .data[[fill_col]],
      sig = !is.na(.data[[sig_col]]) & .data[[sig_col]] < 0.05
    )

  order <- dat %>%
    group_by(Dinuc2) %>%
    summarise(m = mean(fill_val, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(m)) %>%
    pull(Dinuc2)

  dat$Dinuc2 <- factor(dat$Dinuc2, levels = rev(order))
  lim <- max(abs(dat$fill_val), na.rm = TRUE)

  ggplot(dat, aes(Region, Dinuc2, fill = fill_val)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_point(
      data = subset(dat, sig),
      aes(Region, Dinuc2),
      shape = 16, size = 1.6, color = "black"
    ) +
    scale_fill_gradient2(
      low = "#2b8cbe", mid = "white", high = "#d7301f",
      midpoint = 0, limits = c(-lim, lim),
      name = legend_title
    ) +
    scale_x_discrete(labels = function(x) tools::toTitleCase(x)) +
    theme_classic(base_size = 18) +
    theme(
      text = element_text(family = "Arial", colour = "black"),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, colour = "black"),
      axis.text.y = element_text(colour = "black"),
      # plot.title  = element_text(hjust = 0.5, colour = "black"),
      legend.position = "right",
      legend.title = element_text(size = 15, colour = "black"),
      legend.text  = element_text(size = 10, colour = "black"),
      plot.margin = margin(10, 10, 10, 5)
    )
    # + ggtitle(title_text)
}

## ------------------------------------------------------------
## Build panels (each with its own legend)
## ------------------------------------------------------------
p1 <- plot_heatmap(
  all_long,
  "Coefficient_ln_lifespan_c",
  "p.adjusted.fdr_ln_lifespan_c",
  "β1: Lifespan association",
  "Lifespan\nassociation\n(β1)"
)

p2 <- plot_heatmap(
  all_long,
  "Coefficient_ln_Genome_Size_c",
  "p.adjusted.fdr_ln_Genome_Size_c",
  "β2: Genome size association",
  "Genome size\nassociation\n(β2)"
)

p3 <- plot_heatmap(
  all_long,
  "Coefficient_Interaction",
  "p.adjusted.fdr_Interaction",
  "β3: Interaction association",
  "Interaction\nassociation\n(β3)"
)

## ------------------------------------------------------------
## Combine horizontally with ample spacing
## ------------------------------------------------------------
p_combined <- p1 | p2 | p3

ggsave(
  "Fig_heatmap_beta1_beta2_beta3_lnOE_16x6_eachLegend.png",
  p_combined,
  width = 22,
  height = 9,
  dpi = 300
)

cat("Saved: Fig_heatmap_beta1_beta2_beta3_lnOE_16x6_eachLegend.png\n")

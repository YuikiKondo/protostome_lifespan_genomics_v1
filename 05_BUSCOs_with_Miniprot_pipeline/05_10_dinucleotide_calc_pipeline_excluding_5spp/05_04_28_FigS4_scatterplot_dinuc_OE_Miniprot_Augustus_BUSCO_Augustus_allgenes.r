library(ggplot2)
library(dplyr)

# ---- INPUT FILES ----
file_miniprot_busco <- "pgls_results_7_regions_genome10_others16_dinuc_ln_OE_CENTERED_G_CENTERED_interaction_Miniprot_BUSCO_excluding_5spp.csv"
file_augustus_busco <- "pgls_results_7_regions_genome10_others16_dinuc_ln_OE_CENTERED_G_CENTERED_interaction_Augustus_BUSCO_excluding_5spp.csv"
file_augustus_all   <- "pgls_results_7_regions_genome10_others16_dinuc_ln_OE_CENTERED_G_CENTERED_interaction_Augustus_allgenes_excluding_5spp.csv"

# ---- READ ----
miniprot_busco <- read.csv(file_miniprot_busco, check.names = FALSE)
augustus_busco <- read.csv(file_augustus_busco, check.names = FALSE)
augustus_all   <- read.csv(file_augustus_all,   check.names = FALSE)

# ---- KEEP ONLY WHAT WE NEED ----
miniprot_busco <- miniprot_busco %>%
  select(Region, Dinucleotide, Coefficient_lnOE_c) %>%
  rename(beta_miniprot_busco = Coefficient_lnOE_c)

augustus_busco <- augustus_busco %>%
  select(Region, Dinucleotide, Coefficient_lnOE_c) %>%
  rename(beta_augustus_busco = Coefficient_lnOE_c)

augustus_all <- augustus_all %>%
  select(Region, Dinucleotide, Coefficient_lnOE_c) %>%
  rename(beta_augustus_all = Coefficient_lnOE_c)

# ---- MERGE ALL THREE ----
data_merged <- miniprot_busco %>%
  inner_join(augustus_busco, by = c("Region", "Dinucleotide")) %>%
  inner_join(augustus_all,   by = c("Region", "Dinucleotide")) %>%
  filter(Region != "genome")

# ---- SUMS OF ABS BETAS (PER REGION) ----
beta_sums <- data_merged %>%
  group_by(Region) %>%
  summarise(
    Sum_abs_beta_miniprot_busco = sum(abs(beta_miniprot_busco), na.rm = TRUE),
    Sum_abs_beta_augustus_busco = sum(abs(beta_augustus_busco), na.rm = TRUE),
    Sum_abs_beta_augustus_all   = sum(abs(beta_augustus_all),   na.rm = TRUE),
    .groups = "drop"
  )

# ---- PLOTMATH LABELS (atop() for line breaks; no Unicode subscripts) ----
beta_sums_miniprot_vs_aug_busco <- beta_sums %>%
  mutate(
    label_plotmath = paste0(
      "atop(",
      "paste('Sum |', beta[1], '| Miniprot BUSCOs: ', ", round(Sum_abs_beta_miniprot_busco, 2), "),",
      "paste('Sum |', beta[1], '| Augustus BUSCOs: ', ", round(Sum_abs_beta_augustus_busco, 2), ")",
      ")"
    )
  )

beta_sums_all_vs_busco <- beta_sums %>%
  mutate(
    label_plotmath = paste0(
      "atop(",
      "paste('Sum |', beta[1], '| Augustus all genes: ', ", round(Sum_abs_beta_augustus_all, 2), "),",
      "paste('Sum |', beta[1], '| Augustus BUSCOs: ', ", round(Sum_abs_beta_augustus_busco, 2), ")",
      ")"
    )
  )

# ---- GLOBAL AXIS LIMITS (same x/y ranges across all facets and both plots) ----
all_beta_vals <- c(
  data_merged$beta_miniprot_busco,
  data_merged$beta_augustus_busco,
  data_merged$beta_augustus_all
)

lim <- max(abs(all_beta_vals), na.rm = TRUE)
axis_limits <- c(-lim, lim)

# ============================================================
# 1) Augustus ALL genes vs Augustus BUSCOs (beta[1])
# ============================================================
output_file <- "scatterplots_by_region_PGLS_beta_Augustus_all_vs_BUSCOs.png"
png(output_file, width = 1600, height = 1200, res = 150)

ggplot(data_merged, aes(x = beta_augustus_all, y = beta_augustus_busco)) +
  geom_point(size = 2.2, alpha = 0.7, color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(
    ~ Region,
    scales = "fixed",
    labeller = labeller(
      Region = function(x) {
        x <- gsub("([a-zA-Z]+)([0-9]+)", "\\1 \\2", x)  # optional: downstream1 → downstream 1
        tools::toTitleCase(x)
      }
    )
  ) +
  geom_text(
    data = beta_sums_all_vs_busco,
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
  labs(
    x = expression(paste("PGLS coefficient (", beta[1], ") \u2014 Augustus all genes")),
    y = expression(paste("PGLS coefficient (", beta[1], ") \u2014 Augustus BUSCOs"))
  )

dev.off()

# ============================================================
# 2) Miniprot BUSCOs vs Augustus BUSCOs (beta[1])
# ============================================================
output_file <- "scatterplots_by_region_PGLS_beta_Miniprot_vs_Augustus_BUSCOs.png"
png(output_file, width = 1600, height = 1200, res = 150)

ggplot(data_merged, aes(x = beta_miniprot_busco, y = beta_augustus_busco)) +
  geom_point(size = 2.2, alpha = 0.7, color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(
    ~ Region,
    scales = "fixed",
    labeller = labeller(
      Region = function(x) {
        x <- gsub("([a-zA-Z]+)([0-9]+)", "\\1 \\2", x)  # optional: downstream1 → downstream 1
        tools::toTitleCase(x)
      }
    )
  ) +
  geom_text(
    data = beta_sums_miniprot_vs_aug_busco,
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
  labs(
    x = expression(paste("PGLS coefficient (", beta[1], ") \u2014 Miniprot BUSCOs")),
    y = expression(paste("PGLS coefficient (", beta[1], ") \u2014 Augustus BUSCOs"))
  )

dev.off()

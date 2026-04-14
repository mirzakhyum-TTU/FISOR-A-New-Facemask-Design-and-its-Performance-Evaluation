Pareto AUC–ΔP Analysis
# ============================================================#
# INPUT (single file):
#   Penetration_long_with_DeltaP_before_AUC.csv
#
# OUTPUT:
#   AUC_DeltaP_specimen_level.csv
#   AUC_DeltaP_by_design.csv
#   Pareto_AUC_vs_DeltaP.jpg
# ============================================================

# ---- Packages ----
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)

# ---- Read merged dataset (before AUC) ----
dat_long <- read_csv("Penetration_long_with_DeltaP_before_AUC.csv", show_col_types = FALSE) %>%
  mutate(
    SpecimenID      = as.character(SpecimenID),
    Structure       = as.character(Structure),
    Filter          = as.character(Filter),
    ParticleSize_nm = as.numeric(ParticleSize_nm),
    Penetration_pct = as.numeric(Penetration_pct),
    DeltaP_Pa       = as.numeric(DeltaP_Pa)
  ) %>%
  filter(
    !is.na(SpecimenID),
    !is.na(Structure),
    !is.na(Filter),
    !is.na(ParticleSize_nm),
    !is.na(Penetration_pct),
    !is.na(DeltaP_Pa)
  )

# enforce consistent facet & legend ordering
structure_order <- c("Two-layer knit", "Four-layer knit", "Two-layer woven", "Two-layer spunbond")
filter_order    <- c("Activated carbon", "Meltblown", "Glass fiber")

dat_long <- dat_long %>%
  mutate(
    Structure = factor(Structure, levels = structure_order),
    Filter    = factor(Filter, levels = filter_order)
  )

# ---- Compute AUC per specimen (trapezoidal integration over particle size) ----
# AUC units: %·nm (penetration % integrated across size range)
auc_spec <- dat_long %>%
  arrange(SpecimenID, ParticleSize_nm) %>%
  group_by(SpecimenID, Structure, Filter) %>%
  summarise(
    DeltaP_Pa = first(DeltaP_Pa),  # ΔP constant within specimen
    AUC_penetration_pct_nm = {
      x <- ParticleSize_nm
      y <- Penetration_pct
      sum((y[-1] + y[-length(y)]) / 2 * (x[-1] - x[-length(x)]))
    },
    .groups = "drop"
  )

write_csv(auc_spec, "AUC_DeltaP_specimen_level.csv")

# Aggregate to design means (Structure × Filter) ----
design <- auc_spec %>%
  group_by(Structure, Filter) %>%
  summarise(
    mean_AUC = mean(AUC_penetration_pct_nm, na.rm = TRUE),
    mean_DeltaP = mean(DeltaP_Pa, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(design, "AUC_DeltaP_by_design.csv")

# Add labels (S1–S4 and F1–F3) ----
structure_map <- c("Two-layer knit"="S1",
                   "Four-layer knit"="S2",
                   "Two-layer woven"="S3",
                   "Two-layer spunbond"="S4")

filter_map <- c("Activated carbon"="F1",
                "Meltblown"="F2",
                "Glass fiber"="F3")

design <- design %>%
  mutate(
    Design = paste0(structure_map[as.character(Structure)], filter_map[as.character(Filter)])
  )

#  Plot 
shape_map <- c("Activated carbon"=16, "Meltblown"=15, "Glass fiber"=17)

p <- ggplot(design, aes(x = mean_DeltaP, y = mean_AUC)) +
  geom_point(aes(color = Filter, shape = Filter), size = 3.5) +
  geom_text_repel(
    aes(label = Design),
    size = 3.2,
    box.padding = 0.25,
    point.padding = 0.15,
    max.overlaps = Inf,
    min.segment.length = Inf,  # removes leader lines
    segment.color = NA
  ) +
  facet_wrap(~Structure, ncol = 2) +
  scale_shape_manual(values = shape_map) +
  # Requested tick intervals (ticks only; no limits so nothing is clipped)
  scale_x_continuous(breaks = seq(50, 300, by = 50)) +
  scale_y_continuous(breaks = seq(5000, 25000, by = 5000)) +
  labs(
    title = "Filtration–breathability tradeoff (AUC vs ΔP)",
    x = "Pressure drop, ΔP (Pa)",
    y = "Penetration AUC (%·nm)",
    color = "Filter",
    shape = "Filter"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_line(linetype = "dashed", linewidth = 0.4),
    panel.grid.minor = element_blank()
  )

ggsave("Pareto_AUC_vs_DeltaP.jpg", p, width = 11.5, height = 8.0, dpi = 300)

# ============================================================
# END
# ============================================================

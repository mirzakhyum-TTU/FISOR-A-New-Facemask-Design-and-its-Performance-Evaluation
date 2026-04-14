# ==================================================
# Pareto AUC-DeltaP Analysis and qf,o Comparison
# Description:
#   This script computes specimen-level penetration AUC,
#   aggregates design-level Pareto metrics, calculates
#   log10(dp)-weighted mean quality factor (qf,o), and
#   compares Pareto ranks with qf,o ranks.
# ==================================================
 
# ---------------------------
# 1. Load required packages
# ---------------------------
library(dplyr)
library(stringr)
library(tidyr)
 
# ---------------------------
# 2. Define input files
# ---------------------------
pen_file <- "Penetration_long_with_DeltaP_before_AUC.csv"
qf_file  <- "QF.csv"
 
# ---------------------------
# 3. Helper functions
# ---------------------------
 
# Trapezoidal integration for AUC
trapz <- function(x, y) {
  o <- order(x)
  x <- x[o]
  y <- y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
}
 
# Parse specimen IDs of the form SxFy_Rz
parse_id <- function(id) {
  m <- regexec("^S(\\d+)F(\\d+)_R(\\d+)$", id)
  g <- regmatches(id, m)[[1]]
  if (length(g) == 0) return(c(NA, NA, NA))
  c(as.integer(g[2]), as.integer(g[3]), as.integer(g[4]))
}
 
# ---------------------------
# 4. Define structure and filter labels
# ---------------------------
structure_map <- c(
  "Two-layer knit",
  "Four-layer knit",
  "Two-layer woven",
  "Two-layer spunbond"
)
 
filter_map <- c(
  "Activated carbon",
  "Meltblown",
  "Glass fiber"
)
 
# ==================================================
# PART A. Pareto AUC-DeltaP analysis
# ==================================================
 
# ---------------------------
# 5. Read penetration data
# ---------------------------
pen <- read.csv(pen_file, stringsAsFactors = FALSE)
 
pen$ParticleSize_nm <- as.numeric(pen$ParticleSize_nm)
pen$Penetration_pct <- as.numeric(pen$Penetration_pct)
pen$DeltaP_Pa       <- as.numeric(pen$DeltaP_Pa)
pen$SpecimenID      <- trimws(pen$SpecimenID)
 
# ---------------------------
# 6. Compute specimen-level AUC
# ---------------------------
spec_ids <- unique(pen$SpecimenID)
pen_spec <- data.frame()
 
for (id in spec_ids) {
  d <- pen[pen$SpecimenID == id, ]
  d <- d[order(d$ParticleSize_nm), ]
 
  auc <- trapz(d$ParticleSize_nm, d$Penetration_pct)
  dp  <- d$DeltaP_Pa[1]
 
  tmp <- data.frame(
    Specimen_ID    = id,
    DeltaP_Pa      = dp,
    Pen_AUC_pct_nm = auc
  )
 
  pen_spec <- rbind(pen_spec, tmp)
}
 
# ---------------------------
# 7. Add design metadata
# ---------------------------
parsed <- t(sapply(pen_spec$Specimen_ID, parse_id))
 
pen_spec$S <- parsed[, 1]
pen_spec$F <- parsed[, 2]
 
pen_spec$Design    <- paste0("S", pen_spec$S, "F", pen_spec$F)
pen_spec$Structure <- structure_map[pen_spec$S]
pen_spec$Filter    <- filter_map[pen_spec$F]
 
# ---------------------------
# 8. Aggregate to design level
# ---------------------------
pareto <- aggregate(
  cbind(DeltaP_Pa, Pen_AUC_pct_nm) ~ Structure + Filter + Design,
  data = pen_spec,
  FUN = mean,
  na.rm = TRUE
)
 
names(pareto)[names(pareto) == "DeltaP_Pa"]      <- "mean_DeltaP"
names(pareto)[names(pareto) == "Pen_AUC_pct_nm"] <- "mean_PenAUC"
 
# ---------------------------
# 9. Normalize metrics and compute Pareto score
# ---------------------------
pareto$DeltaP_norm <- (pareto$mean_DeltaP - min(pareto$mean_DeltaP)) /
  (max(pareto$mean_DeltaP) - min(pareto$mean_DeltaP))
 
pareto$AUC_norm <- (pareto$mean_PenAUC - min(pareto$mean_PenAUC)) /
  (max(pareto$mean_PenAUC) - min(pareto$mean_PenAUC))
 
pareto$Pareto_score <- sqrt(pareto$DeltaP_norm^2 + pareto$AUC_norm^2)
 
# ==================================================
# PART B. qf,o calculation from QF.csv
# ==================================================
 
# ---------------------------
# 10. Read quality factor data
# ---------------------------
qf <- read.csv(qf_file, stringsAsFactors = FALSE)
qf$Specimen_ID <- trimws(qf$Specimen_ID)
 
# ---------------------------
# 11. Identify qf columns and particle sizes
# ---------------------------
qf_cols <- names(qf)[grepl("^qf_", tolower(names(qf)))]
sizes   <- as.numeric(gsub("\\D+", "", qf_cols))
 
ord     <- order(sizes)
qf_cols <- qf_cols[ord]
sizes   <- sizes[ord]
 
# ---------------------------
# 12. Compute log10(dp) weighting
# ---------------------------
logd  <- log10(sizes)
dlog  <- rep(NA, length(logd))
 
for (i in seq_along(logd)) {
  if (i == 1) {
    dlog[i] <- (logd[i + 1] - logd[i]) / 2
  } else if (i == length(logd)) {
    dlog[i] <- (logd[i] - logd[i - 1]) / 2
  } else {
    dlog[i] <- (logd[i + 1] - logd[i - 1]) / 2
  }
}
 
log_range <- logd[length(logd)] - logd[1]
 
# ---------------------------
# 13. Compute specimen-level qf,o
# ---------------------------
qf_spec <- qf[, c("Specimen_ID", qf_cols)]
 
for (cc in qf_cols) {
  qf_spec[[cc]] <- as.numeric(qf_spec[[cc]])
}
 
qf_spec$qf_o_kPa_inv <- apply(
  qf_spec[, qf_cols],
  1,
  function(v) sum(v * dlog, na.rm = TRUE)
)
 
qf_spec$qf_o_mean_kPa_inv <- qf_spec$qf_o_kPa_inv / log_range
 
# ---------------------------
# 14. Add design metadata to qf,o dataset
# ---------------------------
parsed2 <- t(sapply(qf_spec$Specimen_ID, parse_id))
 
qf_spec$S <- parsed2[, 1]
qf_spec$F <- parsed2[, 2]
 
qf_spec$Design    <- paste0("S", qf_spec$S, "F", qf_spec$F)
qf_spec$Structure <- structure_map[qf_spec$S]
qf_spec$Filter    <- filter_map[qf_spec$F]
 
# ---------------------------
# 15. Compute design-level mean and SD for qf,o
# ---------------------------
qfo_design_mean <- aggregate(
  qf_o_mean_kPa_inv ~ Structure + Filter + Design,
  data = qf_spec,
  mean,
  na.rm = TRUE
)
 
qfo_design_sd <- aggregate(
  qf_o_mean_kPa_inv ~ Structure + Filter + Design,
  data = qf_spec,
  sd,
  na.rm = TRUE
)
 
names(qfo_design_mean)[4] <- "qf_o_mean_kPa_inv_mean"
names(qfo_design_sd)[4]   <- "qf_o_mean_kPa_inv_sd"
 
qfo_design <- merge(
  qfo_design_mean,
  qfo_design_sd,
  by = c("Structure", "Filter", "Design"),
  all = TRUE
)
 
# ==================================================
# PART C. Merge results, rank designs, and export
# ==================================================
 
# ---------------------------
# 16. Merge Pareto and qf,o summaries
# ---------------------------
merged <- merge(
  pareto,
  qfo_design,
  by = c("Structure", "Filter", "Design"),
  all = FALSE
)
 
# ---------------------------
# 17. Rank designs
# ---------------------------
merged$Pareto_rank <- rank(merged$Pareto_score, ties.method = "min")
merged$qfo_rank    <- rank(-merged$qf_o_mean_kPa_inv_mean, ties.method = "min")
 
merged$rank_diff_qfo_minus_pareto <- merged$qfo_rank - merged$Pareto_rank
 
# ---------------------------
# 18. Order final table
# ---------------------------
merged <- merged[order(merged$Pareto_rank), ]
 
# ---------------------------
# 19. Export output files
# ---------------------------
write.csv(merged, "qfo_vs_pareto_comparison.csv", row.names = FALSE)
write.csv(
  merged[order(merged$Design), ],
  "qfo_vs_pareto_dataset_used.csv",
  row.names = FALSE
)
 
# ---------------------------
# 20. Preview final results
# ---------------------------
print(head(merged, 12))

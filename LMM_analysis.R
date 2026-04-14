# --------------------------------------------------
# Linear Mixed-Effects Model (LMM) Analysis
# Study: FISOR Facemask Performance
# Description:
#   Size-dependent particle penetration (30–400 nm)
#   modeled using Structure × Filter × ParticleSize
# Methods:
#   - REML estimation
#   - Type III Satterthwaite tests
# Software:
#   R (version 4.5.2)
# --------------------------------------------------
###  Install Packages 
 install.packages(c("lme4", "lmerTest", "emmeans", "nlme"))
library(lme4)
library(lmerTest)
library(emmeans)
library(nlme)

### Read data
dat <- read.csv("Penetration_LMM_R_long.csv", stringsAsFactors = FALSE)
### Factor coding 
dat$SpecimenID <- factor(dat$SpecimenID)

dat$Structure <- factor(
  dat$Structure,
  levels = c("Two-layer knit", "Four-layer knit", "Two-layer woven", "Two-layer spunbond")
)

dat$Filter <- factor(
  dat$Filter,
  levels = c("Activated carbon", "Meltblown", "Glass fiber")
)

dat$ParticleSize <- factor(dat$ParticleSize_nm, levels = sort(unique(dat$ParticleSize_nm)))

dat$Penetration_pct <- as.numeric(dat$Penetration_pct)

cat("Rows:", nrow(dat), " Unique specimens:", length(unique(dat$SpecimenID)), "\n")

### Create logit-transformed outcome (Sensitivity #1)
eps <- 1e-6
dat$Pen_prop <- pmin(pmax(dat$Penetration_pct/100, eps), 1 - eps)
dat$logitPen <- qlogis(dat$Pen_prop)

### Create Filter×Size variance-group factor for nlme (Sensitivity #2)
dat$Filter_Size <- interaction(dat$Filter, dat$ParticleSize, drop = TRUE)
dat$Filter_Size <- factor(dat$Filter_Size)

############################################################
# PRIMARY MODEL (lmer) - Raw penetration
############################################################
m_full <- lmer(Penetration_pct ~ Structure * Filter * ParticleSize + (1 | SpecimenID),
               data = dat, REML = TRUE)

cat("\n===== PRIMARY MODEL: lmer (raw penetration %) =====\n")
print(summary(m_full))
cat("\n===== ANOVA: m_full =====\n")
print(anova(m_full))

# Residual diagnostics + outlier flags
rp <- resid(m_full, type="pearson")
idx_extreme <- which(abs(rp) > 3)

cat("\nExtreme Pearson residuals (|rp| > 3):", length(idx_extreme), "\n")
if (length(idx_extreme) > 0) {
  print(dat[idx_extreme, c("SpecimenID","Structure","Filter","ParticleSize_nm","Penetration_pct")])
}

# Save diagnostics plots to files
png("Diagnostics_m_full.png", width=1400, height=900)
par(mfrow=c(2,2))
plot(fitted(m_full), resid(m_full), xlab="Fitted", ylab="Residual",
     main="m_full: Residuals vs Fitted"); abline(h=0, lty=2)
qqnorm(resid(m_full), main="m_full: QQ plot"); qqline(resid(m_full))
boxplot(rp ~ dat$ParticleSize, main="m_full: Pearson resid by Size",
        xlab="Particle size", ylab="Pearson resid"); abline(h=c(-3,0,3), lty=c(2,1,2))
boxplot(rp ~ interaction(dat$Structure, dat$Filter), las=2,
        main="m_full: Pearson resid by Structure×Filter",
        xlab="Structure×Filter", ylab="Pearson resid"); abline(h=c(-3,0,3), lty=c(2,1,2))
dev.off()

############################################################
# SENSITIVITY #1 (lmer) - Logit transformed penetration
############################################################
m_logit <- lmer(logitPen ~ Structure * Filter * ParticleSize + (1 | SpecimenID),
                data = dat, REML = TRUE)

cat("\n===== SENSITIVITY #1: lmer (logit penetration) =====\n")
print(summary(m_logit))
cat("\n===== ANOVA: m_logit =====\n")
print(anova(m_logit))

png("Diagnostics_m_logit.png", width=1200, height=600)
par(mfrow=c(1,2))
qqnorm(resid(m_logit), main="m_logit: QQ plot"); qqline(resid(m_logit))
plot(fitted(m_logit), resid(m_logit), xlab="Fitted", ylab="Residual",
     main="m_logit: Residuals vs Fitted"); abline(h=0, lty=2)
dev.off()

############################################################
# SENSITIVITY #2 (nlme) - Heteroscedastic variance model
# BEST MODEL: varIdent by Filter×ParticleSize
############################################################
ctrl <- lmeControl(msMaxIter=200, maxIter=200, niterEM=50, returnObject=TRUE)

# Size-only variance (comparison)
m_var_size <- lme(
  Penetration_pct ~ Structure * Filter * ParticleSize,
  random  = ~ 1 | SpecimenID,
  weights = varIdent(form = ~ 1 | ParticleSize),
  data    = dat,
  method  = "REML",
  control = ctrl
)

# Filter×Size variance (BEST)
m_var_fxsz <- lme(
  Penetration_pct ~ Structure * Filter * ParticleSize,
  random  = ~ 1 | SpecimenID,
  weights = varIdent(form = ~ 1 | Filter_Size),
  data    = dat,
  method  = "REML",
  control = ctrl
)

cat("\n===== SENSITIVITY #2A: nlme (varIdent by Size) =====\n")
print(summary(m_var_size))
cat("\n===== ANOVA: m_var_size =====\n")
print(anova(m_var_size))

cat("\n===== SENSITIVITY #2B: nlme (varIdent by Filter×Size) [BEST] =====\n")
print(summary(m_var_fxsz))
cat("\n===== ANOVA: m_var_fxsz =====\n")
print(anova(m_var_fxsz))

cat("\n===== AIC (REML fits; descriptive) =====\n")
print(AIC(m_var_size, m_var_fxsz))

# Proper comparison for variance structures: ML fits + LRT
m_var_size_ML <- update(m_var_size, method="ML")
m_var_fxsz_ML <- update(m_var_fxsz, method="ML")

cat("\n===== AIC (ML fits) =====\n")
print(AIC(m_var_size_ML, m_var_fxsz_ML))

cat("\n===== LRT: Size variance vs Filter×Size variance (ML) =====\n")
print(anova(m_var_size_ML, m_var_fxsz_ML))

# Diagnostics for nlme (normalized residuals)
png("QQ_m_var_fxsz_normalized.png", width=800, height=800)
qqnorm(resid(m_var_fxsz, type="normalized"), main="m_var_fxsz: QQ (normalized resid)")
qqline(resid(m_var_fxsz, type="normalized"))
dev.off()

############################################################
# Estimated marginal means (EMMs) + Tukey comparisons
# (Use m_full for interpretability on % scale)
############################################################
cat("\n===== EMMs (m_full): Filter within Structure×Size =====\n")
emm_f <- emmeans(m_full, ~ Filter | Structure * ParticleSize)
print(pairs(emm_f, adjust="tukey"))

emm_table <- as.data.frame(emmeans(m_full, ~ Structure * Filter * ParticleSize))
write.csv(emm_table, "EMM_m_full_Structure_Filter_Size.csv", row.names=FALSE)

cat("\nSaved: Diagnostics_m_full.png, Diagnostics_m_logit.png, QQ_m_var_fxsz_normalized.png\n")
cat("Saved: EMM_m_full_Structure_Filter_Size.csv\n")

############################################################
# END
############################################################

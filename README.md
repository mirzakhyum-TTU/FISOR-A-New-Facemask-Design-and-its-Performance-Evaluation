# FISOR-A-New-Facemask-Design-and-its-Performance-Evaluation
R code for analyzing size-dependent particle penetration (30–400 nm) and breathability in FISOR facemasks using mixed-effects modeling and Pareto AUC–ΔP optimization.
# FISOR: A New Facemask Design and Its Performance Evaluation

This repository contains the R code used to analyze size-dependent particle penetration (30–400 nm) and breathability in modular facemask systems (FISOR). The workflow integrates mixed-effects modeling and Pareto AUC–ΔP analysis to evaluate filtration–breathability trade-offs.

---

## Study Overview
The FISOR framework combines multilayer textile supports (S1–S4) with interchangeable filter inserts (F1–F3). A total of 12 structure–filter combinations were evaluated following ASTM F3502-21 protocols.

### Supporting Layers
- S1: Two-layer knit  
- S2: Four-layer knit  
- S3: Two-layer woven  
- S4: Two-layer spunbond  

### Filter Media
- F1: Activated carbon  
- F2: Meltblown  
- F3: Glass fiber  

---

## ⚙️ Methods

### 1. Linear Mixed-Effects Modeling (LMM)
- Response: Particle penetration (%)  
- Fixed effects: Structure × Filter × Particle size  
- Random effects: Specimen ID  
- Estimation: REML  
- Inference: Type III Satterthwaite method  

### 2. Pareto AUC–ΔP Analysis
- AUC: Integrated particle penetration across 30–400 nm (%·nm)  
- ΔP: Pressure drop (Pa)  
- Multi-objective optimization:
  - Minimize penetration (AUC)
  - Minimize breathing resistance (ΔP)
- Ranking based on normalized Euclidean distance from ideal point  

---

## 📂 Repository Structure

- `LMM_analysis.R` → Mixed-effects modeling and statistical analysis  
- `Pareto_AUC_dP_analysis.R` → Pareto ranking and trade-off analysis  
- `data.csv` → Experimental dataset (if included)  
- `README.md` → Project documentation  

---

## ▶️ How to Run

1. Open R (version ≥ 4.0 recommended)  
2. Install required packages:
   install.packages(c("lme4","lmerTest","emmeans","ggplot2","dplyr","tidyr","stringr"))

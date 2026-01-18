# PFAS-Induced-Thyroid-Disruption-Metabolomic-Mediation-Analysis
This repository contains the R scripts used to analyze associations between maternal serum PFAS exposure and thyroid function indices during pregnancy, and to identify metabolomic mediators using MWAS, meet-in-the-middle (MITM) feature selection, and mediation analyses (single-mediator and pathway-level/high-dimensional mediation).

# Study overview
Per- and polyfluoroalkyl substances (PFASs) have been shown to interfere with thyroid function. However, the mechanisms underlying PFAS-induced thyroid disruption in pregnant women remain unclear. This study conducted a comprehensive longitudinal metabolome-wide association study (MWAS) combined with a meet-in-the-middle (MITM) approach to identify key metabolic pathways mediating PFAS-thyroid associations during pregnancy.

- Design: longitudinal pregnancy cohort with repeated measures across two trimesters.
- Exposure: legacy and emerging PFAS
- Outcomes: TSH, FT3, FT4, T3, T4, and T4/T3 ratio.
- Omics: untargeted high-resolution metabolomics
- Main analyses:
  - Linear mixed-effects models for PFAS–thyroid associations
  - Metabolome-wide association study (MWAS) + MITM feature selection
  - LASSO-assisted prioritization
  - Causal mediation analysis (bootstrap)
  - Pathway-level mediation / high-dimensional mediation using `hdmed`
 
# Repository structure
- code/
  - thyroid-metabolites-MWAS_analysis.R
  - sensitivity_analysis-sex.R
  - sensitivity_analysis-imputed.R
  - PFAS-thyroid_analysis.R
  - PFAS-metabolites-MWAS_analysis.R
  - mediation_analysis.R
  -  high-dimensional_mediation.R
  - data_preprocessing.R

- README.md
- requirements.txt
- LICENSE


# Scripts
 - `PFAS-thyroid_analysis.R`  
   Mixed-effects models for PFAS–thyroid outcomes; covariate adjustment; FDR correction.

 - `PFAS-metabolites-MWAS_analysis.R`  
   MWAS: PFAS → metabolomic features (GLM), trimester-specific and/or pooled as designed.

 - `thyroid-metabolites-MWAS_analysis.R`  
   MWAS: metabolomic features → thyroid outcomes (GLM); LASSO-assisted prioritization if used.

 - `mediation_analysis.R`  
   Single-mediator causal mediation analysis (bootstrap) for selected metabolites (e.g., MITM set).

 - `high-dimension_mediation.R`  
   Pathway-level / multi-mediator mediation using `hdmed` (SIS + de-biased LASSO), pathway grouping.

 - `sensitivity_analysis-imputed.R`  
   Sensitivity analysis using multiple imputation for left-censored/LOQ data (distribution-based).

 - `sensitivity_analysis-sex.R`  
   Effect modification / stratified analyses by fetal sex.

# Requirements
 - R >= 4.4.2
 - Install required R packages
    install.packages(c(
      "lme4",           # Mixed-effects models
      "lmerTest",       # Testing in linear mixed models
      "glmnet",         # LASSO regression
      "mediation",      # Mediation analysis
      "hdmed",          # High-dimensional mediation
      "MetaboAnalystR", # Metabolomics pathway analysis
      "ggplot2",        # Visualization
      "pheatmap",       # Heatmaps
      "dplyr",          # Data manipulation
      "tidyr",          # Data tidying
      "reshape2"        # Data reshaping
    ))
  
# Citation
   If you use this code, please cite the associated manuscript.

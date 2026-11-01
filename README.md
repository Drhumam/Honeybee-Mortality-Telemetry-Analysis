# Multi-Sensor Telemetry and Multivariate DTW Clustering for Honeybee Colony Mortality Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

This repository contains the complete R scripts, custom functions, and analysis pipelines developed for the study investigating honeybee colony mortality using multi-sensor hive telemetry, multivariate Dynamic Time Warping (DTW) clustering, machine learning classifiers, and Firth's penalized logistic regression models.

---

## Repository Structure

The workflow is divided into three sequential analytical parts:

```text
├── dataf.csv                             # Raw input telemetry dataset (example / structure)
├── dataf_processed_cleaned.csv           # Imputed and outlier-corrected dataset
├── cluster_stats.csv                     # Summary metrics across DTW clusters (Table 1)
├── 01_preprocessing_dtw.R                # Part 1: Data preprocessing, imputation, and multivariate DTW clustering
├── 02_machine_learning.R                 # Part 2: ML classification pipelines (Random Forest, XGBoost, SVM)
├── 03_shap_firth_models.R                # Part 3: SHAP interpretability analysis and Firth's penalized risk models
├── SHAP_Importance_Combined.png          # High-resolution publication-grade SHAP figure (A & B)
├── Figure_ROC_Comparison_Firth.png       # Publication-grade ROC comparison curve across models
└── LICENSE                               # MIT Open-Source License
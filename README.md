# Ecological assembly dynamics of multi-kingdom gut microbiome predict postoperative infection


This repository contains the custom R scripts used for the analysis presented in the manuscript:

**"Ecological assembly dynamics of multi-kingdom gut microbiome predict postoperative infection"**

This pipeline performs multi-kingdom (Bacteria, Fungi, Archaea) microbiome analysis, including diversity calculation, community state type (CST) clustering, co-occurrence network construction, metabolomics integration, mediation analysis, and machine learning-based prediction of postoperative infections in pediatric liver transplantation (PLT) patients.

## 1. System Requirements

### Software Dependencies
| Software | Version | Required |
|:---------|:--------|:---------|
| R | ≥ 4.0.0 | Yes |
| IQ-TREE (optional) | ≥ 2.3.6 | Only for phylogenetic analysis |

### R Packages
The analysis relies heavily on the following R packages. Please ensure they are installed.

Core data manipulation and visualization
tidyverse ≥ 2.0.0
dplyr ≥ 1.1.0
tidyr ≥ 1.3.0
ggplot2 ≥ 3.4.0
reshape2 ≥ 1.4.4

Ecological and network analysis
vegan ≥ 2.6.0
microeco ≥ 0.15.0
meconetcomp ≥ 0.2.0
igraph ≥ 1.5.0

Machine learning and prediction
caret ≥ 6.0.0
randomForestSRC ≥ 3.2.0
glmnet ≥ 4.1.0
xgboost ≥ 1.7.0
pROC ≥ 1.18.0

Mediation and structural equation modeling
mediation ≥ 4.5.0
lavaan ≥ 0.6.0

Visualization add-ons
ggpubr ≥ 0.6.0
pheatmap ≥ 1.0.12
scatterpie ≥ 0.2.0
UpSetR ≥ 1.4.0



### Operating Systems Tested
| OS | Version | Status |
|:---|:--------|:-------|
| Windows | 10/11 (R ≥ 4.0) |  ✅ Tested |

### Hardware Requirements
- **Minimum**: 4 CPU cores, 16 GB RAM
- **Recommended**: 16+ CPU cores, 64 GB RAM (required for network bootstrap resampling and machine learning hyperparameter tuning on large cohorts)
- **Non-standard hardware**: None required

---

## 2. Installation Guide

### Step 1: Clone Repository
```bash
git clone https://github.com/littlepretty-coder/PLT.git
cd PLT

# Install CRAN packages
packages <- c("tidyverse", "dplyr", "tidyr", "ggplot2", "reshape2", "vegan", "igraph", 
              "caret", "randomForestSRC", "glmnet", "xgboost", "pROC", "mediation", 
              "lavaan", "ggpubr", "pheatmap", "scatterpie", "UpSetR", "RColorBrewer")

install.packages(setdiff(packages, rownames(installed.packages())))

# Install Bioconductor or GitHub packages (if needed)
if (!require("microeco")) install.packages("microeco")
if (!require("meconetcomp")) install.packages("meconetcomp")

Typical Install Time
On a standard desktop computer (4 cores, 16 GB RAM, 100 Mbps internet): ~15-25 minutes.

Main time consumption: Compilation of igraph, xgboost, and tidyverse dependencies.

Typical running Time
On a standard desktop computer (4 cores, 16 GB RAM, 100 Mbps internet): ~15-25 minutes.

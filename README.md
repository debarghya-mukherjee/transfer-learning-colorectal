# Transfer Learning for Colorectal Cancer Survival Prediction

This repository contains R code to replicate experiments for predicting the survival time of rectal cancer patients using gene expression (RNA-Seq TPM) and clinical data from The Cancer Genome Atlas (TCGA). Specifically, the project investigates how data from colon cancer (source domain) can be used to improve prediction performance on rectal cancer patients (target domain), using a fine-tuned ridge regression model based on this paper.

---

## 🔍 Overview

- **Target Task**: Predict log overall survival time (`log(OS.time)`) of rectal cancer patients using gene expression and clinical covariates.
- **Source Data**: Colon cancer data (TCGA-COAD).
- **Target Data**: Rectal cancer data (TCGA-READ).
- **Approach**: Compare a ridge regression model trained only on rectal data ("target only") with a fine-tuned estimator that incorporates predictions from colon cancer data ("transfer learning").

---

## 📁 Directory Structure

```
.
├── transfer_learning_colorectal.R   # Main R script
├── README.md                        # This file
├── .gitignore
└── data/                            # Folder to manually place large downloaded data files
```

---

## 📦 Requirements

Make sure the following R packages are installed:

```r
install.packages(c("glmnet", "dplyr", "tidyr", "Metrics"))
```

---

## 💾 Data Download

Please manually download the following two files and place them in the `data/` directory:

1. `TCGA-COAD.star_tpm.tsv`
2. `TCGA-READ.star_tpm.tsv`

These files are too large to upload directly to GitHub and are available on [Google Drive](https://drive.google.com).

---

## 🚀 Running the Analysis

1. Place the downloaded `.tsv` files into a folder called `data/`.
2. Open `transfer_learning_colorectal.R` in R or RStudio.
3. Update file paths if necessary.
4. Run the script to produce boxplots comparing RMSEs between target-only and fine-tuned estimators.

---

## 📊 Output

A boxplot is generated comparing the prediction performance of the two models (target-only vs. transfer learning) based on RMSE over 500 train/test splits.

---

## 📬 Contact

For questions or suggestions, please contact [Debarghya Mukherjee](mdeb@bu.edu).

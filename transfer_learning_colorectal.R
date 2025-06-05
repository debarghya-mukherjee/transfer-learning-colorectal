# -------------------------------
# Title: Colon-Rectal Cancer Survival Prediction via Transfer Learning
# Author: Debarghya Mukherjee
# Description: Predict log-survival time using gene expression and phenotypic data.
# -------------------------------

# Load required packages
required_packages <- c("glmnet", "dplyr", "tidyr", "Metrics")
lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
})

set.seed(12345)

# -------------------------------
# Load data
# -------------------------------
data_survival <- read.delim("data/TCGA-COAD.survival.tsv")
data_phenotype <- read.delim("data/TCGA-COAD.clinical.tsv")
genome_seq <- read.delim("data/TCGA-COAD.star_tpm.tsv")

data_survival_rectal <- read.delim("data/TCGA-READ.survival.tsv")
genome_rectal_seq <- read.delim("data/TCGA-READ.star_tpm.tsv")
data_phenotype_rectal <- read.delim("data/TCGA-READ.clinical.tsv")

# -------------------------------
# Preprocessing helper
# -------------------------------
id_list_genome = colnames(genome_seq)
id_list_survival = data_survival$sample
id_list_phenotype = data_phenotype$sample

id_list_genome_clean = gsub("\\.", "", id_list_genome)
id_list_survival_clean =  gsub("-", "", id_list_survival)
id_list_phenotype_clean =  gsub("-", "", id_list_phenotype)


# -------------------------------
# Process Colon Data (Source)
# -------------------------------
common_id = intersect(id_list_genome_clean, id_list_survival_clean)

genome_seq_selected = genome_seq[,which(id_list_genome_clean %in% common_id)]
genome_seq_selected = matrix(unlist(genome_seq_selected), nrow(genome_seq_selected), ncol(genome_seq_selected))
genome_seq_selected = t(genome_seq_selected)

data_survival_selected = data_survival$OS.time[which(id_list_survival_clean %in% common_id)]
data_survival_selected = log(as.double(data_survival_selected)) ## log survival

data_phenotype_selected = data_phenotype[which(id_list_phenotype_clean %in% common_id),]


cor_vec = NULL
for(j in 1:ncol(genome_seq_selected))
{
  cor_vec[j] = cor(data_survival_selected, genome_seq_selected[,j])
}

top_idx <- order(abs(cor_vec), decreasing = TRUE, na.last = NA)[1:2500]

selected_cov = c("race.demographic", "age_at_index.demographic", "prior_malignancy.diagnoses")

data_phenotype_selected_reg = data_phenotype_selected[, which(colnames(data_phenotype_selected) %in% selected_cov)]

data_phenotype_selected_reg <- data_phenotype_selected_reg %>%
  mutate(
    is_white = ifelse(race.demographic == "white", 1, 0),
    is_black = ifelse(race.demographic == "black or african american", 1, 0)
  )

data_phenotype_selected_reg$prior_malignancy.diagnoses[data_phenotype_selected_reg$prior_malignancy.diagnoses == "yes"] = 1
data_phenotype_selected_reg$prior_malignancy.diagnoses[data_phenotype_selected_reg$prior_malignancy.diagnoses == "no"] = 0

data_phenotype_selected_reg = data_phenotype_selected_reg %>% select(-c(race.demographic))


# -------------------------------
# Process Rectal Data (Target)
# -------------------------------


id_list_genome_rectal = colnames(genome_rectal_seq)
id_list_survival_rectal = data_survival_rectal$sample
id_list_phenotype_rectal = data_phenotype_rectal$sample

id_list_genome_rectal_clean = gsub("\\.", "", id_list_genome_rectal)
id_list_survival_rectal_clean =  gsub("-", "", id_list_survival_rectal)
id_list_phenotype_rectal_clean =  gsub("-", "", id_list_phenotype_rectal)

common_id_rectal = intersect(id_list_genome_rectal_clean, id_list_survival_rectal_clean)


genome_rectal_seq_selected = genome_rectal_seq[,which(id_list_genome_rectal_clean %in% common_id_rectal)]
genome_rectal_seq_selected = matrix(unlist(genome_rectal_seq_selected), nrow(genome_rectal_seq_selected), ncol(genome_rectal_seq_selected))
genome_rectal_seq_selected = t(genome_rectal_seq_selected)

data_survival_rectal_selected = data_survival_rectal$OS.time[which(id_list_survival_rectal_clean %in% common_id_rectal)]
data_survival_rectal_selected = log(as.double(data_survival_rectal_selected)) ## log survival

cor_vec_rectal = NULL
for(j in 1:ncol(genome_rectal_seq_selected))
{
  cor_vec_rectal[j] = cor(data_survival_rectal_selected, genome_rectal_seq_selected[,j])
}

top_idx_rectal <- order(abs(cor_vec_rectal), decreasing = TRUE, na.last = NA)[1:2500]

common_genes = intersect(top_idx, top_idx_rectal)
print(length(common_genes))

genome_seq_top = genome_seq_selected[,common_genes]
genome_seq_rectal_top = genome_rectal_seq_selected[,common_genes]

data_phenotype_rectal_selected = data_phenotype_rectal[which(id_list_phenotype_rectal_clean %in% common_id_rectal),]


data_phenotype_rectal_selected_reg = data_phenotype_rectal_selected[, which(colnames(data_phenotype_rectal_selected) %in% selected_cov)]

data_phenotype_rectal_selected_reg <- data_phenotype_rectal_selected_reg %>%
  mutate(
    is_white = ifelse(race.demographic == "white", 1, 0),
    is_black = ifelse(race.demographic == "black or african american", 1, 0)
  )

data_phenotype_rectal_selected_reg$prior_malignancy.diagnoses[data_phenotype_rectal_selected_reg$prior_malignancy.diagnoses == "yes"] = 1
data_phenotype_rectal_selected_reg$prior_malignancy.diagnoses[data_phenotype_rectal_selected_reg$prior_malignancy.diagnoses == "no"] = 0

data_phenotype_rectal_selected_reg = data_phenotype_rectal_selected_reg %>% select(-c(race.demographic))


# -------------------------------
# Modeling
# -------------------------------

X_colon = cbind(genome_seq_top, as.matrix(data_phenotype_selected_reg))
Y_colon = data_survival_selected

cv_ridge_source = cv.glmnet(X_colon, Y_colon, alpha = 0, nfolds = 6)


X_rectal = cbind(genome_seq_rectal_top, as.matrix(data_phenotype_rectal_selected_reg))
Y_rectal = data_survival_rectal_selected


Iter = 500
tf_rmse = NULL
target_only_rmse = NULL
n_train <- 100

for (j in 1:Iter) {
  n_train = 100
  train_idx = sample(1:nrow(X_rectal), n_train, replace = F)
  
  X_rectal_train = X_rectal[train_idx,]
  Y_rectal_train = Y_rectal[train_idx]
  
  X_rectal_test = X_rectal[-train_idx,]
  Y_rectal_test = Y_rectal[-train_idx]

  # Target-only Ridge
  
  cv_target_ridge = cv.glmnet(X_rectal_train, Y_rectal_train, alpha = 0, nfolds = 3)
  Y_target_hat = predict(cv_target_ridge, newx = X_rectal_test)
  target_only_rmse[j] = rmse(Y_target_hat, Y_rectal_test)

  # Transfer Learning Ridge
  
  Y_hat_source = predict(cv_ridge_source, newx = X_rectal_train)
  Y_modified = Y_rectal_train - Y_hat_source 
  cv_tf = cv.glmnet(X_rectal_train, Y_modified, alpha = 0, nfolds = 3)
  Y_tf_hat = predict(cv_tf, newx = X_rectal_test)
  Y_tf_hat = Y_tf_hat + predict(cv_ridge_source, newx = X_rectal_test)
  tf_rmse[j] = rmse(Y_tf_hat, Y_rectal_test)
}

# -------------------------------
# Plot
# -------------------------------

pdf("output/rmse_boxplot.pdf", height = 6, width = 6)
boxplot(
  list("Target only" = target_only_rmse, "Transfer" = tf_rmse),
  ylab = "Log survival time RMSE",
  main = "Comparison between target only\nand fine-tuned prediction",
  col = c("#66c2a5", "#fc8d62"),
  border = "black",
  notch = TRUE
)
dev.off()

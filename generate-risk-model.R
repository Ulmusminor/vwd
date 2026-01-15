##############################
# 0. SETUP AND ENVIRONMENT   #
##############################

# Load libraries
library(parallel)
library(qs2)
library(FactoMineR)
library(factoextra)
library(missMDA)
library(mclust)

# Parallel options and random seed
options(rf.cores = detectCores() - 1)
set.seed(1923)

# Working directory
setwd("C:/Users/carlu/Desktop/My_Projects/risk-calculators/von-willebrand")

# Custom functions
source("scripts/preprocess_data.R")

####################################################
# 1. PHENOTYPIC DATA ANALYSIS                      #
####################################################

cat("\n=== PHENOTYPIC DATA ANALYSIS ===\n")

## 1.1 Load data
feno <- read.table("data/Datos_fenotípicos.txt", header = TRUE,
                   sep = "\t", row.names = 1, stringsAsFactors = TRUE, dec = ",")

# Remove "Centro Familia" and "FVIII_B_LC1"
feno$Edat <- as.integer(feno$Edat)
feno <- feno[, -c(1, 13)]

## 1.2 Quality Control
qc_feno <- analizar_dataset(feno)

# Remove column with 95% missing data
feno <- feno[, -12]

## 1.3 Correlation Analysis
corr_feno <- extract_correlations_data(feno)
res_feno <- elim_cor(feno, corr_feno, umbral = 0.98)
feno <- res_feno$evw

## 1.4 FAMD + Clustering
# Impute missing values
feno.imputed <- imputeFAMD(feno, ncp = 3, seed = 1923, maxiter = 1000)$completeObs

# FAMD with 10 and 7 dimensions (compare)
feno.famd <- FAMD(feno.imputed, ncp = 10, graph = TRUE)
fviz_screeplot(feno.famd, addlabels = TRUE)

feno.famd <- FAMD(feno.imputed, ncp = 7, graph = TRUE)
fviz_famd_ind(feno.famd, repel = TRUE, geom = "point")
fviz_famd_var(feno.famd, repel = TRUE)

# Clustering with Mclust
feno.mclust <- Mclust(feno.famd$ind$coord)

cat("\nNúmero óptimo de clusters (fenotípicos):", feno.mclust$G, "\n")
print(table(feno.mclust$classification))
fviz_mclust(feno.mclust, "classification", geom = "point")


####################################################
# 2. GENETIC DATA ANALYSIS                         #
####################################################

cat("\n=== GENETIC DATA ANALYSIS ===\n")

## 2.1 Load data
geno <- read.table("data/Datos_genéticos.txt", header = TRUE,
                   sep = "\t", row.names = 1, dec = ",")

## 2.2 Quality Control
qc_geno <- analizar_dataset(geno)

## 2.3 Correlation Analysis
corr_geno <- extract_correlations_data(geno, umbral = 0.3)
res_geno <- elim_cor(geno, corr_geno, umbral = 0.98)
geno <- res_geno$evw

## 2.4 Remove Frequency Columns
freq_cols <- grep("Freq", names(geno), ignore.case = TRUE, value = TRUE)
geno <- geno[, !(names(geno) %in% freq_cols)]

# QC after filtering
qc_geno_post <- analizar_dataset(geno)

## 2.5 PCA + Clustering
geno.pca <- PCA(geno, ncp = 10, graph = TRUE)
fviz_screeplot(geno.pca, addlabels = TRUE)

geno.pca <- PCA(geno, ncp = 5, graph = TRUE)

geno.mclust <- Mclust(geno.pca$ind$coord)

cat("\nNúmero óptimo de clusters (genéticos):", geno.mclust$G, "\n")
print(table(geno.mclust$classification))
fviz_mclust(geno.mclust, "classification", geom = "point")


####################################################
# 3. SAVE RESULTS (Optional)                       #
####################################################

# You can save QC and clustering results if desired:
# qs_save(list(feno = feno, geno = geno, feno_clusters = feno.mclust, geno_clusters = geno.mclust),
#         "results/vonwillebrand_analysis.qs")

cat("\n=== ANALYSIS COMPLETE ===\n")

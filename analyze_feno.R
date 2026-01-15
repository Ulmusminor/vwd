library(parallel)
library(qs2)
options(rf.cores = detectCores() - 1)
set.seed(1923)
setwd("C:/Users/carlu/Desktop/My_Projects/risk-calculators/von-willebrand")

# Read datasets
feno <- read.table("data/Datos_fenotípicos.txt", header = T,
                   sep = "\t", row.names = 1, stringsAsFactors = T, dec = ",")
# Remove Centro Familia (centro donde son analizados) and FVIII_B_LC1 (todo valor 0)
feno$Edat <- as.integer(feno$Edat)
feno <- feno[,-c(1, 13)]


colnames(feno)

####################################################
# 1. QC                                            #
####################################################
source("scripts/preprocess_data.R")
qc_analysis <- analizar_dataset(feno)

# Remove columna con 95% datos faltantes
feno <- feno[,-12]

####################################################
# 2. Correlation analysis                          #
####################################################
corr_feno <- extract_correlations_data(feno)

# Remove correlated variables
res_feno <- elim_cor(feno, corr_feno, umbral = 0.98)
feno <- res_feno$evw

####################################################
# 3. FAMD Normalization                            #
####################################################
library(FactoMineR)
library(missMDA)
library(factoextra)
feno.imputed <- imputeFAMD(feno, ncp = 3, seed = 1923,
                           maxiter = 1000)$completeObs
feno.famd<- FAMD(feno.imputed, ncp = 10, graph = TRUE)
fviz_screeplot(feno.famd, addlabels = TRUE)
feno.famd<- FAMD(feno.imputed, ncp = 7, graph = TRUE)

fviz_famd_ind(feno.famd, repel = TRUE, geom = "point")  # individuos
fviz_famd_var(feno.famd, repel = TRUE)                  # variables

library(mclust)
res.mclust <- Mclust(feno.famd$ind$coord)

# Número óptimo de clusters
res.mclust$G
# Etiquetas de cluster para cada observación
table(res.mclust$classification)
# Visualizar clusters en el espacio de FAMD
fviz_mclust(res.mclust, "classification", geom = "point")

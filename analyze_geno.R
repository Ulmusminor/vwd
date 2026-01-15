library(parallel)
library(qs2)
library(FactoMineR)
library(mclust)
library(factoextra)

options(rf.cores = detectCores() - 1)
set.seed(1923)

# Read datasets
geno <- read.table("data/Datos_genéticos.txt", header = T,
                   sep = "\t", row.names = 1, dec = ",")

####################################################
# 1. QC                                            #
####################################################
source("scripts/preprocess_data.R")
qc_analysis <- analizar_dataset(geno)

####################################################
# 2. Correlation analysis                          #
####################################################
corr_geno <- extract_correlations_data(geno, umbral = 0.3)

# Remove correlated variables
res_geno <- elim_cor(geno, corr_geno, umbral = 0.98)
geno <- res_geno$evw

####################################################
# 3. Correct FREQ columns                          #
####################################################
freq_cols <- grep("Freq", names(geno), ignore.case = TRUE, value = TRUE)
geno <- geno[, !(names(geno) %in% freq_cols)]

qc_analysis_geno <- analizar_dataset(geno)

####################################################
# 4. FAMD Normalization                            #
####################################################
geno.pca <- PCA(geno, ncp = 10, graph = TRUE)
fviz_screeplot(geno.pca, addlabels = TRUE)
geno.pca <- PCA(geno, ncp = 5, graph = TRUE)

res.mclust <- Mclust(geno.pca$ind$coord)

# Número óptimo de clusters
res.mclust$G
# Etiquetas de cluster para cada observación
table(res.mclust$classification)
# Visualizar clusters en el espacio de FAMD
fviz_mclust(res.mclust, "classification", geom = "point")


# Supongamos:
# - res.mclust$classification son las etiquetas (1..G)
# - hem_score es tu vector numérico con la puntuación hemorrágica
# - opcional: df tiene las covariables (edad, sexo, etc.)

cluster <- factor(res.mclust$classification)
dat <- data.frame(hem_score = feno$Puntuación_hemorrágica, cluster = cluster)

library(ggplot2)
ggplot(dat, aes(cluster, hem_score)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  labs(x = "Cluster (Mclust)", y = "Puntuación hemorrágica")


# ANOVA
fit_aov <- aov(hem_score ~ cluster, data = dat)
summary(fit_aov)

# Efecto (tamaño): eta^2
# install.packages("effectsize")
library(effectsize)
eta_squared(fit_aov)     # η² o partial η²

# Comparaciones post-hoc con ajuste
TukeyHSD(fit_aov)

# install.packages("ordinal")
library(ordinal)
hem_ord <- factor(dat$hem_score, ordered = TRUE)
clm_fit <- clm(hem_ord ~ cluster, data = dat)
summary(clm_fit)



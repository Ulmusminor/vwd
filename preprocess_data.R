analizar_dataset <- function(df) {
  resultados <- data.frame(
    Columna = character(),
    Tipo = character(),
    Detalle = character(),
    Porcentaje_NA = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (col in colnames(df)) {
    columna <- df[[col]]
    tipo <- class(columna)[1]
    porcentaje_na <- mean(is.na(columna)) * 100
    
    if (is.numeric(columna)) {
      detalle <- paste0("Rango: ", 
                        round(min(columna, na.rm = TRUE), 2),
                        " - ", 
                        round(max(columna, na.rm = TRUE), 2))
    } else if (is.factor(columna)) {
      niveles <- paste(levels(columna), collapse = ", ")
      detalle <- paste0("Niveles: ", niveles)
    } else {
      detalle <- "No numérico ni factor"
    }
    
    resultados <- rbind(resultados, data.frame(
      Columna = col,
      Tipo = tipo,
      Detalle = detalle,
      Porcentaje_NA = round(porcentaje_na, 2),
      stringsAsFactors = FALSE
    ))
  }
  
  return(resultados)
}

extract_correlations_data <- function(df, metodo = "spearman", umbral = NULL) {
  # Cargar librerías necesarias
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Por favor instala 'reshape2' con install.packages('reshape2')")
  }
  if (!requireNamespace("corrplot", quietly = TRUE)) {
    stop("Por favor instala 'corrplot' con install.packages('corrplot')")
  }
  
  library(reshape2)
  library(corrplot)
  
  # 1️⃣ Seleccionar solo columnas numéricas
  df_num <- df[, sapply(df, is.numeric)]
  
  if (ncol(df_num) < 2) {
    stop("No hay suficientes columnas numéricas para calcular correlaciones.")
  }
  
  # 2️⃣ Calcular correlación no paramétrica
  cor_mat <- cor(df_num, use = "pairwise.complete.obs", method = metodo)
  
  # 3️⃣ Crear copia para graficar (puede tener NAs si se filtra)
  cor_mat_plot <- cor_mat
  
  # 4️⃣ Convertir a formato largo (tidy)
  cor_df <- melt(cor_mat, varnames = c("Variable1", "Variable2"), value.name = "Correlacion")
  
  # 5️⃣ Eliminar correlaciones consigo mismas
  cor_df <- subset(cor_df, Variable1 != Variable2)
  
  # 6️⃣ Eliminar duplicados simétricos
  cor_df <- cor_df[!duplicated(t(apply(cor_df[, 1:2], 1, sort))), ]
  
  # 7️⃣ Si hay umbral, filtrar correlaciones que lo superen
  if (!is.null(umbral)) {
    cor_df <- subset(cor_df, abs(Correlacion) > umbral)
    
    # También limpiar la matriz para el plot
    cor_mat_plot[abs(cor_mat_plot) <= umbral] <- NA
  }
  
  # 8️⃣ Ordenar por correlación absoluta
  cor_df <- cor_df[order(-abs(cor_df$Correlacion)), ]
  
  # 9️⃣ Graficar matriz de correlaciones
  title_txt <- ifelse(is.null(umbral),
                      paste("Correlaciones (", metodo, ")", sep = ""),
                      paste("Correlaciones (", metodo, ", |ρ| >", umbral, ")", sep = ""))
  
  corrplot(cor_mat_plot,
           method = "color",
           type = "upper",
           tl.cex = 0.6,
           tl.col = "black",
           number.cex = 0.45,
           na.label = " ",
           addCoef.col = "black",
           mar = c(0, 0, 2, 0),
           title = title_txt)
  
  # 10️⃣ Mensaje resumen
  message("Se han encontrado ", nrow(cor_df), " pares de variables correlacionadas.")
  
  # 11️⃣ Devolver el dataframe limpio
  return(cor_df)
}

elim_cor <- function(evw, corr, umbral = 0.98){
  evw <- as.data.frame(evw)
  corr <- subset(corr, !is.na(Correlacion) & abs(Correlacion) > umbral,
                 select = c("Variable1","Variable2","Correlacion"))
  corr$Variable1 <- as.character(corr$Variable1)
  corr$Variable2 <- as.character(corr$Variable2)
  
  drop <- unique(unlist(apply(corr, 1, function(r){
    v1 <- r["Variable1"]; v2 <- r["Variable2"]
    if(v1 %in% names(evw) && v2 %in% names(evw)){
      if(sum(is.na(evw[[v1]])) > sum(is.na(evw[[v2]]))) v1 else v2
    } else NULL
  })))
  
  drop <- intersect(drop, names(evw))
  list(
    evw = evw[, setdiff(names(evw), drop), drop = FALSE],
    eliminadas = drop
  )
}
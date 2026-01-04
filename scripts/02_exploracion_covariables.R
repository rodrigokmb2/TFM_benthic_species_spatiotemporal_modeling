# ======================================================================
# Script 02: Exploración de covariables ambientales y distribución de densidades
# Proyecto: Modelización espaciotemporal de especies bentónicas en ecosistemas costeros de Galicia
# Autor: Rodrigo Alba Salgueiro
# Fecha: 2026-01-01
# ======================================================================

# Librerías
library(tidyverse)   # incluye ggplot2
library(gridExtra)
library(ggcorrplot)


# Cargar datos preparados en 00_setup.R
lances_spp <- readRDS("checkpoints/lances_spp.rds")

# Especies a modelizar
especies_modelo <- c(
  "Polybius henslowii",
  "Pagurus bernhardus",
  "Buglossidium luteum",
  "Pegusa lascaris"
)

cat("Especies a modelizar:\n")
print(especies_modelo)


# Correlación entre covariables continuas
vars_continuas <- lances_spp %>%
  dplyr::select(
    profundidad.x, pctMO, pctArenaGruesa, pctArenaFina,
    pctLodo, Q50phi, So
  ) %>%
  na.omit()   # eliminar filas con NA solo en covariables continuas


cor_matrix   <- cor(vars_continuas, method = "pearson")
cor_pvalues  <- ggcorrplot::cor_pmat(vars_continuas)

# Gráfico
p_cor <- ggcorrplot(
  cor_matrix,
  method   = "circle",
  type     = "lower",
  lab      = TRUE,
  lab_size = 3,
  p.mat    = cor_pvalues,
  sig.level = 0.05
) +
  ggtitle("Matriz de correlación - Covariables continuas")

ggsave(
  "01_correlacion_covariables.png",
  plot   = p_cor,
  width  = 8,
  height = 8,
  dpi    = 300
)

cat("Gráfico guardado: 01_correlacion_covariables.png\n\n")

# Tabla de correlaciones altas (r > 0.7)
high_corr <- cor_matrix
high_corr[lower.tri(high_corr, diag = TRUE)] <- NA

high_corr_table <- as.data.frame(as.table(high_corr)) %>%
  filter(abs(Freq) > 0.7 & !is.na(Freq)) %>%
  arrange(desc(abs(Freq)))

if (nrow(high_corr_table) > 0) {
  cat("Correlaciones > 0.7 detectadas:\n")
  print(high_corr_table)
} else {
  cat("No hay correlaciones > 0.7\n")
}




# Función: gráficos de densidad vs covariables
crear_graficos_densidad <- function(datos, especie_nombre) {
  
  p1 <- ggplot(datos, aes(x = profundidad.x, y = densidad_km2)) +
    geom_point(size = 2, alpha = 0.5, color = "steelblue") +
    geom_smooth(
      method = "loess", se = TRUE,
      color = "red", fill = "red", alpha = 0.2,
      formula = y ~ x
    ) +
    labs(
      title = paste("Densidad vs Profundidad -", especie_nombre),
      x = "Profundidad (m)", 
      y = "Densidad (ind/km²)"
    ) +
    theme_minimal()
  
  p2 <- ggplot(datos, aes(x = zona, y = densidad_km2, fill = zona)) +
    geom_boxplot(alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
    labs(
      title = paste("Densidad por Zona -", especie_nombre),
      x = "Zona", 
      y = "Densidad (ind/km²)"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  p3 <- ggplot(datos, aes(x = sustrato, y = densidad_km2, fill = sustrato)) +
    geom_boxplot(alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
    labs(
      title = paste("Densidad por Sustrato -", especie_nombre),
      x = "Tipo de sustrato", 
      y = "Densidad (ind/km²)"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  p4 <- ggplot(datos, aes(x = ciclo_estacional, y = densidad_km2,
                          fill = ciclo_estacional)) +
    geom_boxplot(alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
    labs(
      title = paste("Densidad por ciclo estacional -", especie_nombre),
      x = "Ciclo estacional", 
      y = "Densidad (ind/km²)"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  p5 <- ggplot(datos, aes(x = Q50phi, y = densidad_km2)) +
    geom_point(size = 2, alpha = 0.5, color = "forestgreen") +
    geom_smooth(
      method = "loess", se = TRUE,
      color = "orange", fill = "yellow", alpha = 0.2,
      formula = y ~ x
    ) +
    labs(
      title = paste("Densidad vs Q50phi -", especie_nombre),
      x = "Q50phi (escala Wentworth)", 
      y = "Densidad (ind/km²)"
    ) +
    theme_minimal()
  
  datos <- datos %>%
    mutate(
      exposicion_lab = dplyr::case_when(
        exposicion == 0 ~ "Protegido",
        exposicion == 1 ~ "Expuesto",
        TRUE            ~ NA_character_
      )
    )
  
  p6 <- ggplot(
    datos,
    aes(
      x    = exposicion_lab,
      y    = densidad_km2,
      fill = exposicion_lab
    )
  ) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.15, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.05, alpha = 0.3, size = 1.5) +
    labs(
      title = paste("Densidad por exposición costera -", especie_nombre),
      x = "Exposición costera", y = "Densidad (ind/km²)"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6)
}


# Función: gráficos de distribución de densidad
crear_graficos_distribucion <- function(datos, especie_nombre) {
  
  p1 <- ggplot(datos, aes(x = densidad_km2)) +
    geom_histogram(bins = 20, fill = "steelblue",
                   color = "black", alpha = 0.7) +
    geom_vline(aes(xintercept = median(densidad_km2, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    labs(
      title = paste("Histograma -", especie_nombre),
      x = "Densidad (ind/km²)", y = "Frecuencia"
    ) +
    theme_minimal()
  
  p2 <- ggplot(datos, aes(x = densidad_km2)) +
    geom_density(fill = "steelblue", alpha = 0.5, color = "black") +
    geom_rug(alpha = 0.3) +
    geom_vline(aes(xintercept = median(densidad_km2, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    labs(
      title = paste("Density plot -", especie_nombre),
      x = "Densidad (ind/km²)", y = "Densidad"
    ) +
    theme_minimal()
  
  datos_no_cero <- datos %>% filter(densidad_km2 > 0)
  
  p3 <- ggplot(datos_no_cero,
               aes(x = log10(densidad_km2 + 0.001))) +
    geom_histogram(bins = 15, fill = "coral",
                   color = "black", alpha = 0.7) +
    geom_vline(
      aes(xintercept = median(log10(densidad_km2 + 0.001),
                              na.rm = TRUE)),
      color = "red", linetype = "dashed", size = 1
    ) +
    labs(
      title = paste("Histograma log10(densidad) -", especie_nombre),
      x = "log10(Densidad + 0.001)", y = "Frecuencia"
    ) +
    theme_minimal()
  
  p4 <- ggplot(datos, aes(sample = densidad_km2)) +
    stat_qq(color = "steelblue", size = 2, alpha = 0.6) +
    stat_qq_line(color = "red", size = 1) +
    labs(title = paste("Q-Q plot -", especie_nombre)) +
    theme_minimal()
  
  list(p1 = p1, p2 = p2, p3 = p3, p4 = p4)
}

# Gráficos de densidad por especie
for (spp in especies_modelo) {
  cat("Procesando:", spp, "\n")
  
  datos_spp <- lances_spp %>%
    filter(especie == spp)
  
  graficos <- crear_graficos_densidad(datos_spp, spp)
  
  grid_plot <- grid.arrange(
    graficos$p1, graficos$p2, graficos$p5,
    graficos$p4, graficos$p3, graficos$p6,
    nrow = 3, ncol = 2
  )
  
  filename <- paste0("02_densidad_", gsub(" ", "_", spp), ".png")
  ggsave(filename, plot = grid_plot, width = 14, height = 10, dpi = 350)
  
  cat("Gráfico guardado:", filename, "\n")
}


# FUNCIÓN: Análisis de ceros y distribución
# Análisis de distribución por especie
analizar_distribucion_gam <- function(datos, especie_nombre) {
  densidad  <- datos$densidad_km2
  n_lances  <- n_distinct(datos$lance)
  n_ceros   <- sum(densidad == 0)
  pct_ceros <- (n_ceros / n_lances) * 100
  
  media    <- mean(densidad, na.rm = TRUE)
  mediana  <- median(densidad, na.rm = TRUE)
  desvest  <- sd(densidad, na.rm = TRUE)
  varianza <- var(densidad, na.rm = TRUE)
  cv       <- (desvest / media) * 100
  
  minimo <- min(densidad, na.rm = TRUE)
  maximo <- max(densidad, na.rm = TRUE)
  
  skewness_val <- mean((densidad - media)^3, na.rm = TRUE) / (desvest^3)
  kurtosis_val <- mean((densidad - media)^4, na.rm = TRUE) / (varianza^2)
  
  data.frame(
    Especie   = especie_nombre,
    N_lances  = n_lances,
    N_ceros   = n_ceros,
    Pct_ceros = round(pct_ceros, 2),
    Media     = round(media, 4),
    Mediana   = round(mediana, 4),
    Desvest   = round(desvest, 4),
    CV_pct    = round(cv, 2),
    Min       = round(minimo, 4),
    Max       = round(maximo, 4),
    Skewness  = round(skewness_val, 4),
    Kurtosis  = round(kurtosis_val, 4)
  )
}

# Tabla 1: estadísticas descriptivas y ceros
estadisticas_gam <- purrr::map_df(
  especies_modelo,
  ~ analizar_distribucion_gam(
    lances_spp %>% filter(especie == .x),
    .x
  )
)

cat("TABLA 1: ESTADÍSTICAS DESCRIPTIVAS Y DIAGNÓSTICO DE CEROS\n")
cat(strrep("-", 120), "\n")
print(estadisticas_gam)
cat(strrep("-", 120), "\n\n")

write.csv(
  estadisticas_gam,
  "05_estadisticas_descriptivas.csv",
  row.names = FALSE
)

cat("Tabla guardada: 05_estadisticas_descriptivas.csv\n\n")


# Tabla 2: diagnóstico para selección de distribución
cat("TABLA 2: DIAGNÓSTICO PARA SELECCIÓN DE DISTRIBUCIÓN\n")
cat(strrep("-", 100), "\n")

diagnostico_distribucion <- estadisticas_gam %>%
  mutate(
    Razon_justificacion = dplyr::case_when(
      Pct_ceros >= 30              ~ "Elevado % ceros → Tweedie",
      CV_pct > 100                 ~ "Alta variabilidad (CV > 100%) → Tweedie",
      Skewness > 1 | Skewness < -1 ~ "Asimetría pronunciada → Tweedie",
      TRUE                         ~ "Distribución aproximadamente normal"
    ),
    Distribucion_sugerida = "Tweedie (GAM)"
  ) %>%
  dplyr::select(
    Especie, Pct_ceros, CV_pct, Skewness,
    Razon_justificacion, Distribucion_sugerida
  )

print(diagnostico_distribucion)
cat(strrep("-", 100), "\n\n")

write.csv(
  diagnostico_distribucion,
  "06_diagnostico_distribucion.csv",
  row.names = FALSE
)

cat("Tabla guardada: 06_diagnostico_distribucion.csv\n\n")


# Gráfico diagnóstico de ceros
diagnostico_data <- estadisticas_gam %>%
  dplyr::select(Especie, Pct_ceros) %>%
  mutate(Especie = gsub(" ", "\n", Especie))

p_diag <- ggplot(diagnostico_data,
                 aes(x = Especie, y = Pct_ceros, fill = Pct_ceros)) +
  geom_col(alpha = 0.7, color = "black", size = 1) +
  geom_hline(yintercept = 30, linetype = "dashed",
             color = "red", linewidth = 1.2) +
  geom_text(
    aes(label = paste0(Pct_ceros, "%")),
    vjust = -0.5, size = 4, fontface = "bold"
  ) +
  scale_fill_gradient(
    low = "lightgreen", high = "darkred", name = "% ceros"
  ) +
  labs(
    title = "Diagnóstico de ceros: justificación distribución Tweedie",
    x = "Especie", y = "Porcentaje de ceros (%)"
  ) +
  theme_minimal()

ggsave(
  "07_diagnostico_ceros.png",
  plot   = p_diag,
  width  = 10,
  height = 6,
  dpi    = 350
)

# Gráficos de distribución por especie
cat("\nGráficos de distribución:\n")

for (spp in especies_modelo) {
  datos_spp <- lances_spp %>%
    filter(especie == spp)
  
  graficos <- crear_graficos_distribucion(datos_spp, spp)
  
  grid_plot <- grid.arrange(
    graficos$p1, graficos$p2, graficos$p3, graficos$p4,
    nrow = 2, ncol = 2
  )
  
  filename <- paste0("04_distribucion_densidad_",
                     gsub(" ", "_", spp), ".png")
  
  ggsave(filename, plot = grid_plot, width = 12, height = 10, dpi = 350)
  cat("Gráfico guardado:", filename, "\n")
}

# Guardar checkpoints
if (!dir.exists("checkpoints")) dir.create("checkpoints")

saveRDS(estadisticas_gam,        "checkpoints/estadisticas_gam.rds")
saveRDS(diagnostico_distribucion,"checkpoints/diagnostico_distribucion.rds")
saveRDS(
  list(
    cor_matrix              = cor_matrix,
    estadisticas_gam        = estadisticas_gam,
    diagnostico_distribucion = diagnostico_distribucion
  ),
  "checkpoints/datos_exploracion_gam.rds"
)

cat("\nExploración de covariables completada.\n")
cat("Tablas guardadas: 05_estadisticas_descriptivas.csv, 06_diagnostico_distribucion.csv\n")

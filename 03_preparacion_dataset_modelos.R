# ======================================================================
# Script 03: Construcción del dataset de modelización (182 lances × 4 especies)
# Proyecto: Modelización espaciotemporal de especies bentónicas en ecosistemas costeros de Galicia
# Autor: Rodrigo Alba Salgueiro
# Fecha: 2026-01-01
# ======================================================================

# Librerías
library(tidyverse)   # incluye dplyr y ggplot2
library(gridExtra)


# Cargar datos preparados en 00_setup.R
lances_spp <- readRDS("checkpoints/lances_spp.rds")

# Estructura básica de lances_spp
cat("Estructura de lances_spp:\n")
cat(sprintf(" - N registros: %d\n", nrow(lances_spp)))
cat(sprintf(" - Lances únicos: %d\n", n_distinct(lances_spp$lance)))
cat(sprintf(" - Especies únicas: %d\n", n_distinct(lances_spp$especie)), "\n\n")

# Extraer variables ambientales por lance (una fila por lance)
lances_variables <- lances_spp %>%
  group_by(lance) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(lance, zona, ciclo_estacional,
    profundidad.x, exposicion, Q50phi, pctMO
  ) %>%
  arrange(lance) %>%
  as_tibble()


cat(sprintf("Total de lances únicos con variables: %d\n\n",
            nrow(lances_variables)))

cat("Verificación de NAs en lances_variables:\n")
print(colSums(is.na(lances_variables)))


# Especies a modelizar (las mismas que en 02 y 04)
especies_modelo <- c(
  "Polybius henslowii",
  "Pagurus bernhardus",
  "Buglossidium luteum",
  "Pegusa lascaris"
)

# Grilla expandida 182 lances × 4 especies
grilla_expandida <- expand_grid(
  lance   = lances_variables$lance,
  especie = especies_modelo
)

cat(sprintf(
  "Grilla expandida: %d lances × %d especies = %d registros\n\n",
  n_distinct(grilla_expandida$lance),
  n_distinct(grilla_expandida$especie),
  nrow(grilla_expandida)
))


# Preparar datos de abundancia desde lances_spp
datos_abundancia <- lances_spp %>%
  dplyr::select(lance, especie, numero_ejemplares, densidad) %>%
  group_by(lance, especie) %>%
  summarise(
    numero_ejemplares = sum(numero_ejemplares, na.rm = TRUE),
    densidad_km2      = mean(densidad, na.rm = TRUE),
    .groups = "drop"
  )

cat(sprintf(
  "Datos de abundancia: %d registros (solo con presencia)\n\n",
  nrow(datos_abundancia)
))

# Unir grilla con abundancia (rellenando ceros)
lancesspp_expandido <- grilla_expandida %>%
  left_join(datos_abundancia, by = c("lance", "especie")) %>%
  mutate(
    numero_ejemplares = replace_na(numero_ejemplares, 0),
    densidad_km2      = replace_na(densidad_km2, 0)
  )

cat(sprintf("Dataset expandido: %d registros\n", nrow(lancesspp_expandido)))
cat(sprintf(
  " - Registros con densidad > 0: %d\n",
  sum(lancesspp_expandido$densidad_km2 > 0)
))
cat(sprintf(
  " - Registros con densidad = 0: %d\n\n",
  sum(lancesspp_expandido$densidad_km2 == 0)
))


# Unir con variables ambientales
lancesspp_completo <- lancesspp_expandido %>%
  left_join(lances_variables, by = "lance")

nas_por_columna <- colSums(is.na(lancesspp_completo))
cat("NAs por columna después de la unión:\n")
print(nas_por_columna)

# Descriptores por especie (antes de filtrar para el modelo)
calcular_descriptores <- function(x) {
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  data.frame(
    n_lances  = length(x),
    presencias = sum(x > 0),
    ceros      = sum(x == 0),
    pct_ceros  = round(mean(x == 0) * 100, 1),
    min        = min(x),
    max        = max(x),
    mediana    = median(x),
    q1         = q1,
    q3         = q3,
    iqr        = q3 - q1,
    q95        = quantile(x, 0.95)
  )
}

descriptores_spp <- lancesspp_expandido %>%
  group_by(especie) %>%
  summarise(across(densidad_km2, calcular_descriptores), .groups = "drop")

if (!dir.exists("checkpoints")) dir.create("checkpoints")

write.csv(descriptores_spp,
          "checkpoints/descriptores_spp.csv",
          row.names = FALSE)
saveRDS(descriptores_spp,
        "checkpoints/descriptores_spp.rds")



# ================================================
# DATASET FINAL PARA MODELIZACIÓN
# ================================================
lancesspp_modelo_completo <- lancesspp_completo %>%
  filter(
    !is.na(profundidad.x),
    !is.na(Q50phi),
    !is.na(zona)
  ) %>%
  mutate(
    densidad_miles   = densidad_km2 / 1000,
    zona             = factor(zona),
    ciclo_estacional = factor(ciclo_estacional),
    exposicion       = factor(exposicion)
  ) %>%
  arrange(especie, lance)

cat("Dataset final para modelización:\n")
cat(sprintf(" - Total registros: %d\n", nrow(lancesspp_modelo_completo)))
cat(sprintf(
  " - Lances únicos: %d\n\n",
  n_distinct(lancesspp_modelo_completo$lance)
))

# Estadísticas por especie
stats_especie <- lancesspp_modelo_completo %>%
  group_by(especie) %>%
  summarise(
    n_obs         = n(),
    n_lances      = n_distinct(lance),
    n_presencias  = sum(densidad_km2 > 0),
    n_ausencias   = sum(densidad_km2 == 0),
    pct_ceros     = round((n_ausencias / n_obs) * 100, 1),
    densidad_media = mean(densidad_km2, na.rm = TRUE),
    densidad_max   = max(densidad_km2, na.rm = TRUE),
    densidad_min   = min(densidad_km2[densidad_km2 > 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_presencias))

print(stats_especie)

write.csv(
  stats_especie,
  "01_estadisticas_especies_182lances.csv",
  row.names = FALSE
)

# Verificación final por especie
especies_lista <- unique(lancesspp_modelo_completo$especie)

cat("\n", strrep("=", 80), "\n")
cat("VERIFICACIÓN FINAL\n")
cat(strrep("=", 80), "\n\n")

for (spp in especies_lista) {
  datos_spp <- lancesspp_modelo_completo %>%
    filter(especie == spp)
  
  cat(spp, ":\n")
  cat(sprintf(
    "  N lances: %d (esperado: 182)\n",
    n_distinct(datos_spp$lance)
  ))
  cat(sprintf(
    "  N registros: %d (esperado: 182)\n",
    nrow(datos_spp)
  ))
  cat(sprintf(
    "  %% ceros: %.1f%%\n",
    (sum(datos_spp$densidad_km2 == 0) / nrow(datos_spp)) * 100
  ))
  cat(sprintf(
    "  Variables ambientales completas: %s\n\n",
    ifelse(sum(is.na(datos_spp$profundidad.x)) == 0, "SÍ", "NO")
  ))
}

# Guardar dataset final (checkpoint para 04)
saveRDS(
  lancesspp_modelo_completo,
  "checkpoints/lancesspp_modelo_completo.rds"
)
write.csv(
  lancesspp_modelo_completo,
  "checkpoints/lancesspp_modelo_completo.csv",
  row.names = FALSE
)

cat("Datos guardados en checkpoints/lancesspp_modelo_completo.[rds,csv]\n")


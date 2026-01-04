# ======================================================================
# Script 00: Preparación de datos y métricas básicas
# Proyecto: Modelización espaciotemporal de especies bentónicas en ecosistemas costeros de Galicia
# Autor: Rodrigo Alba Salgueiro
# Fecha: 2026-01-01
# ======================================================================

# Librerías
library(tidyverse)
library(vegan)
library(BiodiversityR)


# Lectura de datos crudos
datos_spp_raw <- read_csv("Datos/datos_especies.csv")
datos_lances <- read_csv("Datos/datos_lances.csv")
datos_sedimento <- read_delim("Datos/datos_sedimento.csv", delim = ";")

# Renombrar columnas de estación
      # estación numérica >>> (estacion)
      # estación categórica -fría/cálida- >>> ciclo_estacional
datos_lances <- datos_lances %>%
  rename(
    estacion = estacion...2,
    ciclo_estacional = estacion...13
  )

# Limpieza antes de unir:
# Identificar lances válidos (datos_lances)
lances_validos <- unique(datos_lances$lance)

# Filtrar datos de especies a lances válidos
datos_spp <- datos_spp_raw %>%
  filter(lance %in% lances_validos)

# Convertir a enteros
datos_lances <- datos_lances %>%
  mutate(
    estacion = as.integer(estacion),
    lance    = as.integer(lance)
  )

# Agregar especies por lance, y estación
datos_spp <- datos_spp %>%
  group_by(estacion, lance, cod, comercial, nombre_comun,
           nombre_Galicia, FAO_code, especie) %>%
  summarise(
    numero_ejemplares = sum(numero_ejemplares),
    .groups = "drop"
  )

# Unir datos de especies con variables de lances
lances_spp <- datos_spp %>%
  inner_join(
    datos_lances %>%
      dplyr::select(
        lance, estacion, mes, ciclo_estacional, zona, sustrato,
        profundidad, area_arrastrada, exposicion, lat_firmes, lon_firmes
      ),
    by = c("lance", "estacion")
  )

# Unir con datos de sedimento
lances_spp <- lances_spp %>%
  left_join(
    datos_sedimento %>% dplyr::select(-zona),
    by = "estacion"
  )

# Calcular densidad
lances_spp <- lances_spp %>%
  mutate(
    densidad     = numero_ejemplares / area_arrastrada,
    densidad_km2 = (numero_ejemplares / area_arrastrada) * 1000
  )

# Comprobación:
# n_distinct(datos_spp_raw$lance)   # Lances brutos: 209 (en datos_spp_raw)
# n_distinct(datos_spp$lance)       # Lances válidos: 182 (en datos_spp)


# Exploración básica
n_lances <- n_distinct(datos_spp$lance)
n_spp    <- n_distinct(datos_spp$especie)

cat("Número total de lances:", n_lances, "\n")
cat("Número total de especies:", n_spp, "\n")

# Comprobar diferencia  de lances entre datos_spp VS lances_spp
lances_lances_spp <- sort(unique(lances_spp$lance))
lances_datos_spp  <- sort(unique(datos_spp$lance))

lances_extra <- setdiff(lances_datos_spp, lances_lances_spp)
length(lances_extra)
lances_extra

# Número de spp distintas por estación numérica
spp_por_estacion <- datos_spp %>%
  group_by(estacion) %>%
  summarise(n_spp = n_distinct(especie), .groups = "drop")

print(spp_por_estacion)

# Número de spp distintas por ciclo estacional (fría/cálida)
spp_por_ciclo <- datos_spp %>%
  inner_join(
    dplyr::select(datos_lances, estacion, lance, ciclo_estacional),
    by = c("estacion", "lance")) %>%
  group_by(ciclo_estacional) %>%
  summarise(n_spp = n_distinct(especie), .groups = "drop")

print(spp_por_ciclo)


# Frecuencia de aparición de spp por lance
frecuencia_spp <- datos_spp %>%
  group_by(especie) %>%
  summarise(
    lances_con_presenza = n_distinct(lance),
    .groups = "drop"
  ) %>%
  mutate(
    proporcion_lances = lances_con_presenza / n_lances
  ) %>%
  filter(proporcion_lances >= 0.2) %>%  # ≥ 20 % de los lances
  arrange(desc(proporcion_lances), desc(lances_con_presenza))

print(frecuencia_spp)


# Matriz de abundancias por estación (filas = estaciones, columnas = especies)
abund_by_station <- datos_spp %>%
  group_by(estacion, especie) %>%
  summarise(n = sum(numero_ejemplares), .groups = "drop") %>%
  pivot_wider(
    names_from = especie,
    values_from = n,
    values_fill = 0
  ) %>%
  column_to_rownames("estacion")

# 1.b Matriz de abundancias por zona (filas = zonas, columnas = especies)
abund_by_zona <- lances_spp %>%
  group_by(zona, especie) %>%
  summarise(n = sum(numero_ejemplares), .groups = "drop") %>%
  pivot_wider(
    names_from = especie,
    values_from = n,
    values_fill = 0
  ) %>%
  column_to_rownames("zona")

# Abundancia total por especie
abundancia_total <- datos_spp %>%
  group_by(especie) %>%
  summarise(
    total = sum(numero_ejemplares),
    .groups = "drop"
  ) %>%
  arrange(desc(total)) %>%
  mutate(rank = row_number())

# Total de lances
total_lances <- n_distinct(lances_spp$lance)
cat(sprintf("Total de lances en el dataset: %d\n\n", total_lances))

# Métricas por especie
metricas_spp <- lances_spp %>%
  group_by(especie) %>%
  summarise(
    lances_presencia = n_distinct(lance),
    abundancia_total = sum(numero_ejemplares, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    abundancia_relativa_pct =
      (abundancia_total / sum(abundancia_total)) * 100,
    ocurrencia_pct =
      (lances_presencia / total_lances) * 100
  )

# Tabla enriquecida
abundancia_total_enriquecida <- abundancia_total %>%
  left_join(
    metricas_spp %>%
      dplyr::select(especie,
             lances_presencia,
             abundancia_relativa_pct,
             ocurrencia_pct),
    by = "especie"
  ) %>%
  arrange(rank)

# Guardar CSV
write.csv(
  abundancia_total_enriquecida,
  "abundancia_total_enriquecida.csv",
  row.names = FALSE
)

cat("\nTabla guardada: abundancia_total_enriquecida.csv\n")
cat(sprintf("Total de especies: %d\n", nrow(abundancia_total_enriquecida)))
cat(sprintf("Total de individuos: %.0f\n",
            sum(abundancia_total_enriquecida$total)))
cat(sprintf("Total de lances: %d\n\n", total_lances))

# Checkpoints
if (!dir.exists("checkpoints")) dir.create("checkpoints")

# Dataset principal
saveRDS(lances_spp, "checkpoints/lances_spp.rds")

# Matrices y métricas usadas en 01, 02, 03
saveRDS(
  list(
    abund_by_station              = abund_by_station,
    abund_by_zona                 = abund_by_zona,
    abundancia_total              = abundancia_total,
    metricas_spp                  = metricas_spp,
    abundancia_total_enriquecida  = abundancia_total_enriquecida,
    spp_por_estacion              = spp_por_estacion,
    spp_por_ciclo                 = spp_por_ciclo,
    frecuencia_spp                = frecuencia_spp,
    total_lances                  = total_lances
  ),
  "checkpoints/matrices_y_metricas.rds"
)
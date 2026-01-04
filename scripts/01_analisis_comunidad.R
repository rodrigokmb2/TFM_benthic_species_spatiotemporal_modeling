# ======================================================================
# Script 01: Análisis de comunidad (diversidad, ANOVA, nMDS y SIMPER)
# Proyecto: Modelización espaciotemporal de especies bentónicas en ecosistemas costeros de Galicia
# Autor: Rodrigo Alba Salgueiro
# Fecha: 2026-01-01
# ======================================================================

# Librerías
library(tidyverse)
library(vegan)

# Cargar datos del checkpoint 00
lances_spp <- readRDS("checkpoints/lances_spp.rds")
matrices   <- readRDS("checkpoints/matrices_y_metricas.rds")

abund_by_station             <- matrices$abund_by_station
abund_by_zona                <- matrices$abund_by_zona
abundancia_total             <- matrices$abundancia_total
metricas_spp                 <- matrices$metricas_spp
abundancia_total_enriquecida <- matrices$abundancia_total_enriquecida
spp_por_estacion             <- matrices$spp_por_estacion
spp_por_ciclo                <- matrices$spp_por_ciclo
frecuencia_spp               <- matrices$frecuencia_spp
total_lances                 <- matrices$total_lances

# Función para índices de diversidad
calc_diversity <- function(mat) {
  data.frame(
    shannon = diversity(mat, index = "shannon"),
    simpson = diversity(mat, index = "simpson"),
    riqueza = rowSums(mat > 0),
    pielou  = diversity(mat, index = "shannon") / log(rowSums(mat > 0))
  )
}


# Índices por estación y zona
indices_estacion <- calc_diversity(as.matrix(abund_by_station)) %>%
  tibble::rownames_to_column("estacion")

indices_zona <- calc_diversity(as.matrix(abund_by_zona)) %>%
  tibble::rownames_to_column("zona")

print(indices_estacion)
print(indices_zona)


# Curva de acumulación de especies por lance
mat_com_lance <- lances_spp %>%
  group_by(lance, especie) %>%
  summarise(n = sum(numero_ejemplares), .groups = "drop") %>%
  pivot_wider(names_from = especie, values_from = n, values_fill = 0)

mat_lances <- mat_com_lance %>%
  column_to_rownames("lance") %>%
  as.matrix()

set.seed(123)
sac_lances <- specaccum(mat_lances, method = "random")

plot(
  sac_lances,
  ci.type = "poly",
  col = rgb(70/255, 130/255, 180/255, 0.4),
  border = NA,
  xlab = "Número de lances",
  ylab = "Número acumulado de especies",
  main = "Curva de acumulación de especies por lance"
)
lines(sac_lances, col = "steelblue", lwd = 2)
points(sac_lances$sites, sac_lances$richness, pch = 21, bg = "steelblue")


# Normalidad de índices (para ANOVA)
shapiro.test(indices_estacion$shannon)
shapiro.test(indices_estacion$simpson)
shapiro.test(indices_estacion$riqueza)
shapiro.test(indices_estacion$pielou)

# Nota: todos los índices cumplieron normalidad (p > 0.05),
# por lo que se aplicó ANOVA paramétrico por zona.

# Crear el dataframe para relacionar las estaciones y las zonas
estacion_zona <- lances_spp %>%
  dplyr::select(estacion, zona) %>%
  distinct()

indices_estacion <- indices_estacion %>%
  dplyr::mutate(estacion = as.character(estacion)) %>%
  left_join(
    estacion_zona %>% dplyr::mutate(estacion = as.character(estacion)),
    by = "estacion"
  )

# ANOVA para cada índice, comparando las 3 zonas
anova_shannon_zona <- aov(shannon ~ zona, data = indices_estacion)
anova_simpson_zona <- aov(simpson ~ zona, data = indices_estacion)
anova_riqueza_zona <- aov(riqueza ~ zona, data = indices_estacion)
anova_pielou_zona  <- aov(pielou ~ zona,  data = indices_estacion)

summary(anova_shannon_zona)
summary(anova_simpson_zona)
summary(anova_riqueza_zona)
summary(anova_pielou_zona)



###############
# ANÁLISIS nMDS
###############

# Análisis nMDS por estación y zona
# Matriz de comunidades: filas = estaciones, columnas = especies
mat_comunidades <- lances_spp %>%
  group_by(estacion, especie) %>%
  summarise(n = sum(numero_ejemplares), .groups = "drop") %>%
  pivot_wider(
    names_from = especie,
    values_from = n,
    values_fill = 0)

meta_estaciones <- lances_spp %>%
  dplyr::select(estacion, zona) %>%
  distinct()

mat_nmds <- mat_comunidades %>%
  column_to_rownames("estacion") %>%
  as.matrix()

dist_bc <- vegdist(mat_nmds, method = "bray")

set.seed(123)
nmds_res <- metaMDS(dist_bc, k = 2, trymax = 100)

cat("Stress nMDS:", nmds_res$stress, "\n")



# Grafico nMDS, color por zona y con ordihull
zonas <- meta_estaciones$zona[match(rownames(mat_nmds), meta_estaciones$estacion)]
zonas <- factor(zonas)
col_zonas <- as.numeric(zonas)

plot(
  nmds_res, type = "n",
  main = "Ordenación nMDS de la comunidad bentónica (Bray–Curtis)"
)

ordihull(
  nmds_res,
  groups = zonas,
  draw   = "lines",
  col    = col_zonas,
  lwd    = 2.5
)

points(
  nmds_res,
  display = "sites",
  pch     = 21,
  bg      = col_zonas,
  col     = "lightgrey",
  cex     = 2
)

legend(
  "topright",
  legend = levels(zonas),
  pt.bg  = 1:length(levels(zonas)),
  pch    = 21,
  bty    = "n"
)

mtext(
  paste("Stress =", round(nmds_res$stress, 3)),
  side = 3, line = 0.5, cex = 1.0
)



########
# SIMPER
########

simper_zonas <- simper(mat_nmds, group = zonas, permutations = 999)
simper_sum <- summary(simper_zonas)

# Guardar cada contraste en un data frame separado
simper_ac_or <- as.data.frame(simper_sum[[1]])
simper_ac_ab <- as.data.frame(simper_sum[[2]])
simper_or_ab <- as.data.frame(simper_sum[[3]])

head(simper_ac_ab)
head(simper_ac_or)
head(simper_or_ab)

simper_ac_ab_sig <- subset(simper_ac_ab, p < 0.05)
simper_ac_or_sig <- subset(simper_ac_or, p < 0.05)
simper_or_ab_sig <- subset(simper_or_ab, p < 0.05)

# Guardar resultados de comunidad
if (!dir.exists("checkpoints")) dir.create("checkpoints")

saveRDS(
  list(
    indices_estacion   = indices_estacion,
    indices_zona       = indices_zona,
    anova_shannon_zona = anova_shannon_zona,
    anova_simpson_zona = anova_simpson_zona,
    anova_riqueza_zona = anova_riqueza_zona,
    anova_pielou_zona  = anova_pielou_zona,
    nmds_res           = nmds_res,
    mat_nmds           = mat_nmds,
    meta_estaciones    = meta_estaciones,
    simper_zonas       = simper_zonas,
    simper_ac_or       = simper_ac_or,
    simper_ac_ab       = simper_ac_ab,
    simper_or_ab       = simper_or_ab,
    simper_ac_or_sig   = simper_ac_or_sig,
    simper_ac_ab_sig   = simper_ac_ab_sig,
    simper_or_ab_sig   = simper_or_ab_sig,
    sac_lances         = sac_lances,
    mat_lances         = mat_lances
  ),
  "checkpoints/resultados_comunidad.rds"
)

cat("Análisis de comunidad completado.\n")
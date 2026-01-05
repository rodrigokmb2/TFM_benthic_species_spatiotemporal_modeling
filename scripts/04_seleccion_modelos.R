# ======================================================================
# Script 04: Selección de modelos y gráficos QQ y residuos
# Proyecto: Modelización espaciotemporal de especies bentónicas en ecosistemas costeros de Galicia
# Autor: Rodrigo Alba Salgueiro
# Fecha: 2026-01-01
# ======================================================================

# Librerías
library(tidyverse)
library(mgcv)
library(MASS)
library(gridExtra)
library(DHARMa)
library(mgcViz)

# Cargar datos preparados en 03_preparar_182_lances.R
lancesspp_modelo_completo <- readRDS("checkpoints/lancesspp_modelo_completo.rds")

# Filtrar datos completos para modelización
lancesspp_modelo <- lancesspp_modelo_completo %>%
  filter(
    !is.na(densidad_km2),
    !is.na(Q50phi),
    !is.na(profundidad.x),
    !is.na(zona),
    !is.na(pctMO)
  )

especies_modelo <- c(
  "Polybius henslowii",
  "Pagurus bernhardus",
  "Buglossidium luteum",
  "Pegusa lascaris"
)

# ------------------------------------------
# Función para ajustar y comparar modelos
# ------------------------------------------
ajustar_modelos <- function(datos, especie_nombre) {
  
  cat("\n", strrep("=", 70), "\n")
  cat("Especie:", especie_nombre, "\n")
  cat(strrep("=", 70), "\n\n")
  
  datos_spp <- datos %>%
    dplyr::filter(especie == especie_nombre)
  
  n_obs     <- nrow(datos_spp)
  n_lances  <- n_distinct(datos_spp$lance)
  pct_ceros <- (sum(datos_spp$densidad_km2 == 0) / n_obs) * 100
  
  cat(sprintf("N obs: %d | N lances: %d | %% ceros: %.1f%%\n\n",
              n_obs, n_lances, pct_ceros))
  
  # Modelo 1: Negative Binomial
  modelo_nb <- MASS::glm.nb(
    densidad_km2 ~ profundidad.x + Q50phi + pctMO +
      zona + ciclo_estacional + exposicion,
    data = datos_spp
  )
  
  # Modelo 2: Tweedie GAM con potencia fija (p = 1.5)
  modelo_tw_fijo <- mgcv::gam(
    densidad_km2 ~ s(profundidad.x, k = 4) +
      s(Q50phi,        k = 4) +
      s(pctMO,         k = 4) +
      zona + ciclo_estacional + exposicion,
    family = mgcv::tw(theta = 1.5),
    data   = datos_spp,
    method = "REML"
  )
  
  # Modelo 3: Tweedie GAM con potencia estimada
  modelo_tw_estimado <- mgcv::gam(
    densidad_km2 ~ s(profundidad.x, k = 4) +
      s(Q50phi,        k = 4) +
      s(pctMO,         k = 4) +
      zona + ciclo_estacional + exposicion,
    family = mgcv::tw(),
    data   = datos_spp,
    method = "REML"
  )
  
  # Modelo 4: Poisson
  modelo_poisson <- glm(
    densidad_km2 ~ profundidad.x + Q50phi + pctMO +
      zona + ciclo_estacional + exposicion,
    family = poisson(link = "log"),
    data   = datos_spp
  )
  
  extraer_criterios <- function(mod) {
    aic <- AIC(mod)
    bic <- BIC(mod)
    dev <- if (!is.null(mod$deviance) && !is.null(mod$null.deviance)) {
      (1 - mod$deviance / mod$null.deviance) * 100
    } else {
      NA
    }
    c(AIC = aic, BIC = bic, Deviance = dev)
  }
  
  crit_nb          <- extraer_criterios(modelo_nb)
  crit_tw_fijo     <- extraer_criterios(modelo_tw_fijo)
  crit_tw_estimado <- extraer_criterios(modelo_tw_estimado)
  crit_pois        <- extraer_criterios(modelo_poisson)
  
  comparacion <- tibble(
    Modelo   = c("NB", "Tweedie (p = 1.5)", "Tweedie (p estimado)", "Poisson"),
    AIC      = c(crit_nb["AIC"],          crit_tw_fijo["AIC"],
                 crit_tw_estimado["AIC"], crit_pois["AIC"]),
    BIC      = c(crit_nb["BIC"],          crit_tw_fijo["BIC"],
                 crit_tw_estimado["BIC"], crit_pois["BIC"]),
    Deviance = c(crit_nb["Deviance"],     crit_tw_fijo["Deviance"],
                 crit_tw_estimado["Deviance"], crit_pois["Deviance"])
  )
  
  print(comparacion)
  cat("\n")
  
  list(
    especie          = especie_nombre,
    n_obs            = n_obs,
    n_lances         = n_lances,
    pct_ceros        = pct_ceros,
    comparacion      = comparacion,
    modelo_nb        = modelo_nb,
    modelo_tw_fijo   = modelo_tw_fijo,
    modelo_tw_estimado = modelo_tw_estimado,
    modelo_poisson   = modelo_poisson
  )
}

# ---------------------------------------
# Ajustar modelos para todas las especies
# ---------------------------------------

resultados_modelos <- vector("list", length(especies_modelo))
names(resultados_modelos) <- especies_modelo

for (spp in especies_modelo) {
  resultados_modelos[[spp]] <- ajustar_modelos(lancesspp_modelo, spp)
}

# --------------------------------
# Función de diagnóstico DHARMa # GAM Tweedie p estimado
# --------------------------------
plot_dharma_diagnostics <- function(modelo, especie_nombre,
                                    sufijo = "",
                                    nsim = 1000) {
  
  sim <- DHARMa::simulateResiduals(
    fittedModel = modelo,
    n = nsim
  )
  
  filename_diag <- paste0(
    "07_diagnostico_DHARMa_",
    gsub(" ", "_", especie_nombre),
    sufijo,
    ".png"
  )
  
  png(filename_diag, width = 2400, height = 1200, res = 300)
  
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
  
  # 1) QQ + uniformidad
  DHARMa::plotQQunif(sim, main = "QQ plot (uniformidad)")
  abline(0, 1, col = "red", lwd = 2)
  
  # 2) Residual vs predicted (interno de DHARMa)
  DHARMa::plotResiduals(sim, main = "Residual vs predicted")
  
  dev.off()
  
  # Tests SOLO en consola (como debe ser)
  cat("\n", especie_nombre, "\n")
  print(DHARMa::testUniformity(sim))
  print(DHARMa::testDispersion(sim))
  print(DHARMa::testOutliers(sim))
  
  cat("Diagnóstico DHARMa guardado:", filename_diag, "\n")
}

# ----------------------------------------------------------------------
# Generar diagnósticos DHARMa para el modelo con p estimado
# ----------------------------------------------------------------------

for (spp in especies_modelo) {
  plot_dharma_diagnostics(
    modelo = resultados_modelos[[spp]]$modelo_tw_estimado,
    especie_nombre = spp,
    sufijo = "_tweedie_p_estimado",
    nsim = 1000
  )
}

# Efectos de variables
par(mfrow = c(2, 2))
plot(resultados_modelos[["Polybius henslowii"]]$modelo_tw_estimado,
     pages = 1, residuals = TRUE)
par(mfrow = c(2, 2))

plot(resultados_modelos[["Pagurus bernhardus"]]$modelo_tw_estimado,
     pages = 1, residuals = TRUE)
par(mfrow = c(2, 2))

plot(resultados_modelos[["Buglossidium luteum"]]$modelo_tw_estimado,
     pages = 1, residuals = TRUE)
par(mfrow = c(2, 2))

plot(resultados_modelos[["Pegusa lascaris"]]$modelo_tw_estimado,
     pages = 1, residuals = TRUE)



# -----------------------------------------------------
# Diagnóstico básico GAM (Q-Q + residuos vs ajustados)
# -----------------------------------------------------
plot_gam_basic_diagnostics <- function(modelo, especie_nombre,
                                         sufijo = "") {
  
  # Residuos tipo deviance
  res <- residuals(modelo, type = "deviance")
  
  # Valores ajustados en escala de la respuesta
  fitted_vals <- fitted(modelo)
  
  filename_diag <- paste0(
    "07_diagnostico_basico_GAM_",
    gsub(" ", "_", especie_nombre),
    sufijo,
    ".png"
  )
  
  png(filename_diag, width = 2400, height = 1200, res = 300)
  
  par(mfrow = c(1, 2),
      mar = c(5, 5, 4, 2),
      cex = 0.95)
  
  # ---------------------------
  # Q-Q plot de residuos
  # ---------------------------
  stats::qqnorm(
    res,
    main = "Q-Q plot de residuos (deviance)",
    pch  = 16,
    cex  = 0.7
  )
  stats::qqline(res, col = "red", lwd = 2)
  
  # --------------------------------
  # Residuos vs valores ajustados
  # --------------------------------
  plot(
    fitted_vals, res,
    xlab = "Valores ajustados",
    ylab = "Residuos (deviance)",
    main = "Residuos vs valores ajustados",
    pch  = 16,
    cex  = 0.7
  )
  abline(h = 0, col = "red", lwd = 2)
  
  dev.off()
  
  cat("Diagnóstico básico GAM guardado:", filename_diag, "\n")
}

# -----------------------------------------------------
# Bucle diagnósticos básicos para GAM Tweedie p estimado
# -----------------------------------------------------
for (spp in especies_modelo) {
  
  plot_gam_basic_diagnostics(
    modelo = resultados_modelos[[spp]]$modelo_tw_estimado,
    especie_nombre = spp,
    sufijo = "_tweedie_p_estimado"
  )
}


# -----------------------------------------------------
# Tabla resumen de modelos y guardado de checkpoints
# -----------------------------------------------------

# Función auxiliar para extraer criterios de un modelo
extraer_criterios <- function(mod) {
  aic <- AIC(mod)
  bic <- BIC(mod)
  dev <- if (!is.null(mod$deviance) && !is.null(mod$null.deviance)) {
    (1 - mod$deviance / mod$null.deviance) * 100
  } else {
    NA
  }
  c(AIC = aic, BIC = bic, Deviance = dev)
}

tabla_final <- do.call(
  rbind,
  lapply(resultados_modelos, function(x) {
    
    crit_nb          <- extraer_criterios(x$modelo_nb)
    crit_tw_fijo     <- extraer_criterios(x$modelo_tw_fijo)
    crit_tw_estimado <- extraer_criterios(x$modelo_tw_estimado)
    crit_pois        <- extraer_criterios(x$modelo_poisson)
    
    comparacion <- tibble::tibble(
      Modelo   = c("NB", "Tweedie (p = 1.5)", "Tweedie (p estimado)", "Poisson"),
      AIC      = c(crit_nb["AIC"],          crit_tw_fijo["AIC"],
                   crit_tw_estimado["AIC"], crit_pois["AIC"]),
      BIC      = c(crit_nb["BIC"],          crit_tw_fijo["BIC"],
                   crit_tw_estimado["BIC"], crit_pois["BIC"]),
      Deviance = c(crit_nb["Deviance"],     crit_tw_fijo["Deviance"],
                   crit_tw_estimado["Deviance"], crit_pois["Deviance"])
    ) %>%
      dplyr::arrange(AIC)
    
    mejor      <- comparacion$Modelo[1]
    aic_mejor  <- comparacion$AIC[1]
    bic_mejor  <- comparacion$BIC[1]
    dev_mejor  <- comparacion$Deviance[1]
    
    data.frame(
      Especie       = x$especie,
      N_obs         = x$n_obs,
      N_lances      = x$n_lances,
      Pct_ceros     = round(x$pct_ceros, 1),
      Modelo_optimo = mejor,
      AIC           = round(aic_mejor, 2),
      BIC           = round(bic_mejor, 2),
      Deviance_pct  = round(dev_mejor, 2),
      stringsAsFactors = FALSE
    )
  })
)

cat("\n", strrep("=", 80), "\n")
cat("RESUMEN COMPARATIVO - CRITERIOS DE INFORMACIÓN\n")
cat(strrep("=", 80), "\n\n")

print(tabla_final)

write.csv(
  tabla_final,
  "07_resumen_modelos_comparativo.csv",
  row.names = FALSE
)

cat("\nTabla guardada: 07_resumen_modelos_comparativo.csv\n\n")

# Guardar checkpoints necesarios para el script 05
if (!dir.exists("checkpoints")) dir.create("checkpoints")

saveRDS(resultados_modelos,
        "checkpoints/resultados_modelos.rds")
saveRDS(lancesspp_modelo,
        "checkpoints/lancesspp_modelo_filtrado_para_modelos.rds")

cat("Checkpoints guardados.\n")


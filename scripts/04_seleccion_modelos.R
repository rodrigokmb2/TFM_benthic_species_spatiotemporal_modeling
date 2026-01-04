# ======================================================================
# Script 04: Modelización de densidad con GAM-Tweedie y modelos alternativos
# Proyecto: Comunidad bentónica y modelos GAM-Tweedie (Golfo Ártabro)
# Autor: Rodrigo Alba Salgueiro
# Fecha: 2026-01-01
# ======================================================================

# Librerías
library(tidyverse)
library(mgcv)
library(MASS)
library(gridExtra)

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

# Función de diagnóstico gráfico para modelos GAM
plot_gam_diagnostics <- function(modelo, especie_nombre) {
  if (!inherits(modelo, "gam")) {
    cat("El modelo para", especie_nombre, "no es un GAM; se omite diagnóstico.\n")
    return(invisible(NULL))
  }
  
  # Extraer residuos y valores ajustados
  res   <- residuals(modelo, type = "deviance")
  ajust <- fitted(modelo)
  
  df_diag <- tibble(
    ajustados = ajust,
    residuos  = res
  )
  
  p_qq <- ggplot(df_diag, aes(sample = residuos)) +
    stat_qq(color = "steelblue", size = 1.5, alpha = 0.7) +
    stat_qq_line(color = "red", linewidth = 1) +
    labs(
      title = paste("Q-Q residuos devianza -", especie_nombre),
      x     = "Cuantiles teóricos",
      y     = "Residuos de devianza"
    ) +
    theme_minimal()
  
  p_resid <- ggplot(df_diag, aes(x = ajustados, y = residuos)) +
    geom_point(alpha = 0.5, color = "steelblue", size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = paste("Residuos vs valores ajustados -", especie_nombre),
      x     = "Valores ajustados",
      y     = "Residuos de devianza"
    ) +
    theme_minimal()
  
  grid_diag <- gridExtra::grid.arrange(p_qq, p_resid, nrow = 1)
  
  filename_diag <- paste0(
    "07_diagnostico_residuos_",
    gsub(" ", "_", especie_nombre),
    ".png"
  )
  
  ggsave(
    filename_diag,
    plot   = grid_diag,
    width  = 12,
    height = 5,
    dpi    = 300
  )
  
  cat("Diagnóstico de residuos guardado:", filename_diag, "\n")
}


# Función para ajustar los cuatro modelos por especie
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
  ) %>%
    arrange(AIC) %>%
    mutate(
      DeltaAIC = AIC - min(AIC),
      Rank     = dplyr::row_number()
    )
  
  cat("Comparación de modelos:\n")
  print(comparacion)
  cat("\n")
  
  mejor_modelo <- comparacion$Modelo[1]
  aic_mejor    <- comparacion$AIC[1]
  bic_mejor    <- comparacion$BIC[1]
  dev_mejor    <- comparacion$Deviance[1]
  
  cat(sprintf(
    "Modelo seleccionado: %s (AIC = %.2f, BIC = %.2f, Deviance = %.2f%%)\n\n",
    mejor_modelo, aic_mejor, bic_mejor, dev_mejor
  ))
  
  modelo_final <- switch(
    mejor_modelo,
    "NB"                   = modelo_nb,
    "Tweedie (p = 1.5)"    = modelo_tw_fijo,
    "Tweedie (p estimado)" = modelo_tw_estimado,
    "Poisson"              = modelo_poisson
  )
  
  cat(strrep("-", 70), "\n")
  cat("Resumen del modelo elegido\n")
  cat(strrep("-", 70), "\n\n")
  
  if (inherits(modelo_final, "gam")) {
    suma <- summary(modelo_final)
    print(suma)
    
    # Coeficientes paramétricos
    if (!is.null(suma$p.table) && nrow(suma$p.table) > 0) {
      tabla_param <- as.data.frame(suma$p.table)
      tabla_param$Term <- rownames(tabla_param)
      colnames(tabla_param) <- c(
        "Estimate", "Std.Error", "z.value", "p.value", "Term"
      )
      tabla_param <- tabla_param %>%
        dplyr::select(Term, Estimate, Std.Error, z.value, p.value)
      
      filename_param <- paste0(
        "07_parametros_", gsub(" ", "_", especie_nombre), ".csv"
      )
      write.csv(tabla_param, filename_param, row.names = FALSE)
      cat("Tabla de parámetros:", filename_param, "\n")
    }
    
    # Smooth terms
    if (!is.null(suma$s.table) && nrow(suma$s.table) > 0) {
      tabla_smooth <- as.data.frame(suma$s.table)
      tabla_smooth$Term <- rownames(tabla_smooth)
      colnames(tabla_smooth) <- c(
        "edf", "Ref.df", "F", "p.value", "Term"
      )
      tabla_smooth <- tabla_smooth %>%
        dplyr::select(Term, edf, Ref.df, F, p.value)
      
      filename_smooth <- paste0(
        "07_smooth_", gsub(" ", "_", especie_nombre), ".csv"
      )
      write.csv(tabla_smooth, filename_smooth, row.names = FALSE)
      cat("Tabla de smooth terms:", filename_smooth, "\n")
    }
  } else {
    print(coef(summary(modelo_final)))
  }

  # Rangos y niveles para predicciones
  rango_prof <- range(datos_spp$profundidad.x, na.rm = TRUE)
  rango_q50  <- range(datos_spp$Q50phi, na.rm = TRUE)
  
  lvl_zona  <- levels(datos_spp$zona)
  lvl_ciclo <- levels(datos_spp$ciclo_estacional)
  lvl_expo  <- levels(datos_spp$exposicion)
  
  # Profundidad
  pred_prof <- expand.grid(
    profundidad.x    = seq(rango_prof[1], rango_prof[2], length.out = 50),
    Q50phi           = median(datos_spp$Q50phi, na.rm = TRUE),
    pctMO            = median(datos_spp$pctMO, na.rm = TRUE),
    zona             = lvl_zona[1],
    ciclo_estacional = lvl_ciclo[1],
    exposicion       = lvl_expo[1]
  )
  pred_prof$pred <- predict(modelo_final, newdata = pred_prof, type = "response")
  
  p1 <- ggplot(pred_prof, aes(x = profundidad.x, y = pred)) +
    geom_line(color = "steelblue", size = 1) +
    geom_point(
      data = datos_spp,
      aes(x = profundidad.x, y = densidad_km2),
      alpha = 0.3, size = 1
    ) +
    labs(
      title = paste(especie_nombre, "- Profundidad"),
      x = "Profundidad (m)", y = "Densidad (ind/km²)"
    ) +
    theme_minimal()
  
  # Q50phi
  pred_q50 <- expand.grid(
    profundidad.x    = median(datos_spp$profundidad.x, na.rm = TRUE),
    Q50phi           = seq(rango_q50[1], rango_q50[2], length.out = 50),
    pctMO            = median(datos_spp$pctMO, na.rm = TRUE),
    zona             = lvl_zona[1],
    ciclo_estacional = lvl_ciclo[1],
    exposicion       = lvl_expo[1]
  )
  pred_q50$pred <- predict(modelo_final, newdata = pred_q50, type = "response")
  
  p2 <- ggplot(pred_q50, aes(x = Q50phi, y = pred)) +
    geom_line(color = "darkgreen", size = 1) +
    geom_point(
      data = datos_spp,
      aes(x = Q50phi, y = densidad_km2),
      alpha = 0.3, size = 1
    ) +
    labs(
      title = paste(especie_nombre, "- Q50phi"),
      x = "Q50phi", y = "Densidad (ind/km²)"
    ) +
    theme_minimal()
  
  # Zona
  pred_zona <- expand.grid(
    profundidad.x    = median(datos_spp$profundidad.x, na.rm = TRUE),
    Q50phi           = median(datos_spp$Q50phi, na.rm = TRUE),
    pctMO            = median(datos_spp$pctMO, na.rm = TRUE),
    zona             = lvl_zona,
    ciclo_estacional = lvl_ciclo[1],
    exposicion       = lvl_expo[1]
  )
  pred_zona$pred <- predict(modelo_final, newdata = pred_zona, type = "response")
  
  p3 <- ggplot(pred_zona, aes(x = zona, y = pred, fill = zona)) +
    geom_col(alpha = 0.7) +
    labs(
      title = paste(especie_nombre, "- Zona"),
      x = "Zona", y = "Densidad (ind/km²)"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Ciclo estacional
  pred_ciclo <- expand.grid(
    profundidad.x    = median(datos_spp$profundidad.x, na.rm = TRUE),
    Q50phi           = median(datos_spp$Q50phi, na.rm = TRUE),
    pctMO            = median(datos_spp$pctMO, na.rm = TRUE),
    zona             = lvl_zona[1],
    ciclo_estacional = lvl_ciclo,
    exposicion       = lvl_expo[1]
  )
  pred_ciclo$pred <- predict(modelo_final, newdata = pred_ciclo, type = "response")
  
  p4 <- ggplot(
    pred_ciclo,
    aes(x = ciclo_estacional, y = pred, fill = ciclo_estacional)
  ) +
    geom_col(alpha = 0.7) +
    labs(
      title = paste(especie_nombre, "- Ciclo estacional"),
      x = "Ciclo estacional", y = "Densidad (ind/km²)"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  grid_efectos <- grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
  filename_fig <- paste0("07_efectos_", gsub(" ", "_", especie_nombre), ".png")
  
  ggsave(
    filename_fig,
    plot   = grid_efectos,
    width  = 14,
    height = 10,
    dpi    = 300
  )
  
  cat("Figura de efectos:", filename_fig, "\n")
  
  # Gráficos de diagnóstico del modelo final (si es GAM)
  plot_gam_diagnostics(modelo_final, especie_nombre)
  
  list(
    especie        = especie_nombre,
    n_obs          = n_obs,
    n_lances       = n_lances,
    pct_ceros      = pct_ceros,
    comparacion    = comparacion,
    mejor_modelo   = mejor_modelo,
    aic_mejor      = aic_mejor,
    bic_mejor      = bic_mejor,
    deviance_mejor = dev_mejor,
    modelo_objeto  = modelo_final
  )
}

# ----------------------------------------------------------------------
# Ajustar modelos para todas las especies
# ----------------------------------------------------------------------

resultados_tweedie <- vector("list", length(especies_modelo))
names(resultados_tweedie) <- especies_modelo

for (spp in especies_modelo) {
  resultados_tweedie[[spp]] <- ajustar_modelos(lancesspp_modelo, spp)
}

# Resumen comparativo global
cat("\n", strrep("=", 80), "\n")
cat("RESUMEN COMPARATIVO - CRITERIOS DE INFORMACIÓN\n")
cat(strrep("=", 80), "\n\n")

tabla_final <- do.call(
  rbind,
  lapply(resultados_tweedie, function(x) {
    data.frame(
      Especie       = x$especie,
      N_obs         = x$n_obs,
      N_lances      = x$n_lances,
      Pct_ceros     = round(x$pct_ceros, 1),
      Modelo_optimo = x$mejor_modelo,
      AIC           = round(x$aic_mejor, 2),
      BIC           = round(x$bic_mejor, 2),
      Deviance_pct  = round(x$deviance_mejor, 2),
      stringsAsFactors = FALSE
    )
  })
)

print(tabla_final)

write.csv(
  tabla_final,
  "07_resumen_modelos_comparativo.csv",
  row.names = FALSE
)

cat("\nTabla guardada: 07_resumen_modelos_comparativo.csv\n\n")

# Guardar checkpoints
if (!dir.exists("checkpoints")) dir.create("checkpoints")

saveRDS(resultados_tweedie,
        "checkpoints/resultados_tweedie.rds")
saveRDS(lancesspp_modelo,
        "checkpoints/lancesspp_modelo_filtrado_para_modelos.rds")

cat("Checkpoints guardados.\n")

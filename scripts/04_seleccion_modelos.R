# ======================================================================
# Script 04: Selección de modelos y gráficos Q-Q y residuos
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
  dplyr::filter(
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

# Funciones helpers
extraer_criterios <- function(mod) {
  aic <- AIC(mod)
  bic <- BIC(mod)
  dev <- if (!is.null(mod$deviance) && !is.null(mod$null.deviance)) {
    (1 - mod$deviance / mod$null.deviance) * 100
  } else {
    NA_real_
  }
  c(AIC = aic, BIC = bic, Deviance = dev)
}

safe_name <- function(x) gsub(" ", "_", x)


# ===============================
# AJUSTE Y COMPARACIÓN DE MODELOS
# ===============================
# Función para ajustar y comparar modelos
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

# Bucle de ajuse de modelos
resultados_modelos <- vector("list", length(especies_modelo))
names(resultados_modelos) <- especies_modelo

for (spp in especies_modelo) {
  resultados_modelos[[spp]] <- ajustar_modelos(lancesspp_modelo, spp)
}

# =========================================
# Diagnóstico DHARMa GAM Tweedie p estimado
# =========================================
set.seed(123)

for (spp in especies_modelo) {

  mod <- resultados_modelos[[spp]]$modelo_tw_estimado
  if (is.null(mod)) next

  cat("\n", strrep("-", 60), "\n")
  cat("Diagnóstico DHARMa –", spp, "\n")
  cat(strrep("-", 60), "\n")

  sim <- DHARMa::simulateResiduals(
    fittedModel = mod,
    n = 1000
  )

  # Plots
  DHARMa::plotQQunif(
    sim,
    main = paste("DHARMa – QQ uniformidad:", spp)
  )

  DHARMa::plotResiduals(sim)  # sin main
  mtext(paste("DHARMa – Residuos vs predicho:", spp),
        side = 3, line = 0.5, cex = 0.85, adj = 0)

  # Tests solo texto
  print(DHARMa::testUniformity(sim))
  print(DHARMa::testDispersion(sim))
  print(DHARMa::testOutliers(sim))
}

# ====================================================
# Diagnóstico básico GAM (Q-Q + residuos vs ajustados)
# ====================================================
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
  
  png(filename_diag, width = 2400, height = 1200, res = 350)
  
  par(mfrow = c(1, 2),
      mar = c(5, 5, 4, 2),
      cex = 0.95)
  
  # =============================
  # Q-Q plot de residuos
  # =============================
  stats::qqnorm(
    res,
    main = "Q-Q plot de residuos (deviance)",
    pch  = 16,
    cex  = 0.7
  )
  stats::qqline(res, col = "red", lwd = 2)
  
  # =============================
  # Residuos vs valores ajustados
  # =============================
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

# Bucle Plot diagnósticos básicos para GAM Tweedie p estimado
for (spp in especies_modelo) {
  
  plot_gam_basic_diagnostics(
    modelo = resultados_modelos[[spp]]$modelo_tw_estimado,
    especie_nombre = spp,
    sufijo = "_tweedie_p_estimado"
  )
}

# ============================================
# Plot efectos parciales GAM Tweedie p estimado
# ============================================
vars_cont <- c("profundidad.x", "Q50phi", "pctMO")
vars_fact <- c("zona", "ciclo_estacional", "exposicion")

# Colores sobrios (académicos)
col_line  <- "#2C7FB8"  # azul
col_ribbon <- "#2C7FB8" # mismo color (alpha en ribbon)
col_point <- "#D95F0E"  # naranja sobrio

mode_level <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  tb <- table(x)
  names(tb)[which.max(tb)]
}

predict_ci_resp <- function(modelo, newdata) {
  pr <- predict(modelo, newdata = newdata, type = "link", se.fit = TRUE)
  fit <- as.numeric(pr$fit)
  se  <- as.numeric(pr$se.fit)
  
  tibble(
    fit_link = fit,
    lwr_link = fit - 1.96 * se,
    upr_link = fit + 1.96 * se,
    fit_resp = modelo$family$linkinv(fit),
    lwr_resp = modelo$family$linkinv(fit - 1.96 * se),
    upr_resp = modelo$family$linkinv(fit + 1.96 * se)
  )
}

build_ref_row <- function(datos_spp) {
  
  ref_prof <- median(datos_spp$profundidad.x, na.rm = TRUE)
  ref_q50  <- median(datos_spp$Q50phi, na.rm = TRUE)
  ref_mo   <- median(datos_spp$pctMO, na.rm = TRUE)
  
  ref_zona  <- mode_level(datos_spp$zona)
  ref_ciclo <- mode_level(datos_spp$ciclo_estacional)
  ref_expo  <- mode_level(datos_spp$exposicion)
  
  tibble(
    profundidad.x    = ref_prof,
    Q50phi           = ref_q50,
    pctMO            = ref_mo,
    zona             = if (is.factor(datos_spp$zona)) factor(ref_zona, levels = levels(datos_spp$zona)) else factor(as.character(ref_zona)),
    ciclo_estacional = if (is.factor(datos_spp$ciclo_estacional)) factor(ref_ciclo, levels = levels(datos_spp$ciclo_estacional)) else factor(as.character(ref_ciclo)),
    exposicion       = if (is.factor(datos_spp$exposicion)) factor(ref_expo, levels = levels(datos_spp$exposicion)) else factor(as.character(ref_expo))
  )
}

build_pred_cont_1spp <- function(modelo, datos_spp, varname) {
  
  x <- datos_spp[[varname]]
  x_seq <- seq(
    from = as.numeric(quantile(x, 0.02, na.rm = TRUE)),
    to   = as.numeric(quantile(x, 0.98, na.rm = TRUE)),
    length.out = 200
  )
  
  nd0 <- build_ref_row(datos_spp)
  nd  <- nd0[rep(1, length(x_seq)), ]
  nd[[varname]] <- x_seq
  
  pr <- predict_ci_resp(modelo, nd)
  bind_cols(tibble(x = x_seq), pr)
}

plot_cont_1spp <- function(df, varname) {
  ggplot(df, aes(x = x, y = fit_resp)) +
    geom_ribbon(aes(ymin = lwr_resp, ymax = upr_resp),
                fill = col_ribbon, alpha = 0.20) +
    geom_line(color = col_line, linewidth = 1) +
    labs(
      title = varname,
      x = varname,
      y = "Densidad predicha (ind/km²)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.title = element_text(size = 9),
      axis.text  = element_text(size = 8)
    )
}

build_pred_fact_1spp <- function(modelo, datos_spp, varname) {
  
  nd0 <- build_ref_row(datos_spp)
  
  x <- datos_spp[[varname]]
  lvls <- if (is.factor(x)) levels(x) else sort(unique(as.character(x)))
  lvls <- lvls[!is.na(lvls)]
  
  nd <- nd0[rep(1, length(lvls)), ]
  nd[[varname]] <- factor(lvls, levels = lvls)
  
  pr <- predict_ci_resp(modelo, nd)
  tibble(level = nd[[varname]]) %>% bind_cols(pr)
}

plot_fact_1spp <- function(df, varname) {
  ggplot(df, aes(x = level, y = fit_resp)) +
    geom_pointrange(aes(ymin = lwr_resp, ymax = upr_resp),
                    color = col_point) +
    labs(
      title = varname,
      x = varname,
      y = "Densidad predicha (ind/km²)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.title = element_text(size = 9),
      axis.text  = element_text(size = 8),
      axis.text.x = element_text(angle = 25, hjust = 1)
    )
}

# ===============================
# Generar figuras por especie 2x3
# ===============================
for (spp in especies_modelo) {
  
  modelo <- resultados_modelos[[spp]]$modelo_tw_estimado
  datos_spp <- lancesspp_modelo %>% dplyr::filter(especie == spp)
  
  if (is.null(modelo) || nrow(datos_spp) == 0) next
  
  p_cont <- lapply(vars_cont, function(v) {
    dfv <- build_pred_cont_1spp(modelo, datos_spp, v)
    plot_cont_1spp(dfv, v)
  })
  
  p_fact <- lapply(vars_fact, function(v) {
    dfv <- build_pred_fact_1spp(modelo, datos_spp, v)
    plot_fact_1spp(dfv, v)
  })
  
  titulo <- grid::textGrob(
    paste0("Efectos parciales – GAM Tweedie (p estimado): ", spp),
    gp = grid::gpar(fontface = "bold", fontsize = 12)
  )
  
  fig <- gridExtra::arrangeGrob(
    grobs = c(p_cont, p_fact),
    nrow = 2, ncol = 3,
    top = titulo
  )
  
  out_file <- paste0("08_efectos_parciales_GAMTweedie_", gsub(" ", "_", spp), ".png")
  ggsave(out_file, fig, width = 13, height = 7.5, dpi = 350)
  
  cat("✓ Guardado:", out_file, "\n")
}

cat("\n Figuras generadas en el directorio raíz.\n")

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


# ======================================================================
# Script 05: Ajuste del parámetro k y de la complejidad de smooths en GAM-Tweedie
# Proyecto: Modelización espaciotemporal de especies bentónicas en ecosistemas costeros de Galicia
# Autor: Rodrigo Alba Salgueiro
# Fecha: 2026-01-01
# ======================================================================

# Librerías
library(tidyverse)
library(mgcv)

# Cargar datos preparados
if (file.exists("checkpoints/lancesspp_modelo_filtrado_para_modelos.rds")) {
  lancesspp_modelo <- readRDS("checkpoints/lancesspp_modelo_filtrado_para_modelos.rds")
} else if (file.exists("checkpoints/lancesspp_modelo_completo.rds")) {
  lancesspp_modelo <- readRDS("checkpoints/lancesspp_modelo_completo.rds") %>%
    filter(
      !is.na(densidad_km2),
      !is.na(Q50phi),
      !is.na(profundidad.x),
      !is.na(zona),
      !is.na(pctMO)
    )
} else {
  stop("Falta checkpoints/lancesspp_modelo_completo.rds (script 03).")
}

lancesspp_modelo <- lancesspp_modelo %>%
  mutate(
    zona             = if (is.factor(zona)) zona else factor(zona),
    ciclo_estacional = if (is.factor(ciclo_estacional)) ciclo_estacional else factor(ciclo_estacional),
    exposicion       = if (is.factor(exposicion)) exposicion else factor(exposicion)
  )

especies_modelo <- c(
  "Polybius henslowii",
  "Pagurus bernhardus",
  "Buglossidium luteum",
  "Pegusa lascaris"
)


# A) Tuning de k en MODELO 3 (GAM Tweedie p estimado)
k_grid <- c(3, 5)

ajustar_modelo3_k <- function(datos_spp, k) {
  mgcv::gam(
    densidad_km2 ~
      s(profundidad.x, k = k) +
      s(Q50phi,        k = k) +
      s(pctMO,         k = k) +
      zona + ciclo_estacional + exposicion,
    family = mgcv::tw(),
    data   = datos_spp,
    method = "REML"
  )
}


extraer_info_modelo <- function(mod) {
  
  dev_pct <- (1 - mod$deviance / mod$null.deviance) * 100
  theta   <- mod$family$getTheta(TRUE)
  
  suma <- summary(mod)
  
  edf_prof <- edf_q50 <- edf_mo <- NA_real_
  if (!is.null(suma$s.table) && nrow(suma$s.table) > 0) {
    st <- as.data.frame(suma$s.table)
    st$Term <- rownames(st)
    
    pick_edf <- function(pat) {
      ix <- which(grepl(pat, st$Term, fixed = TRUE))
      if (length(ix) == 0) return(NA_real_)
      as.numeric(st$edf[ix[1]])
    }
    
    edf_prof <- pick_edf("profundidad.x")
    edf_q50  <- pick_edf("Q50phi")
    edf_mo   <- pick_edf("pctMO")
  }
  
  k_index_prof <- kp_prof <- k_index_q50 <- kp_q50 <- k_index_mo <- kp_mo <- NA_real_
  
  kchk <- mgcv::k.check(mod)
  kc <- as.data.frame(kchk)
    kc$Term <- rownames(kc)
    
    grab_k <- function(pat) {
      ix <- which(grepl(pat, kc$Term, fixed = TRUE))
      if (length(ix) == 0) return(c(NA_real_, NA_real_))
      kcol <- grep("k", names(kc), ignore.case = TRUE, value = TRUE)
      pcol <- grep("p", names(kc), ignore.case = TRUE, value = TRUE)
      k_val <- if (length(kcol) > 0) as.numeric(kc[ix[1], kcol[1]]) else NA_real_
      p_val <- if (length(pcol) > 0) as.numeric(kc[ix[1], pcol[1]]) else NA_real_
      c(k_val, p_val)
    }
    
    kp <- grab_k("profundidad.x")
    kq <- grab_k("Q50phi")
    km <- grab_k("pctMO")
    
    k_index_prof <- kp[1];  kp_prof <- kp[2]
    k_index_q50  <- kq[1];  kp_q50 <- kq[2]
    k_index_mo   <- km[1];  kp_mo  <- km[2]

  
  tibble(
    AIC          = AIC(mod),
    BIC          = BIC(mod),
    Deviance_pct = dev_pct,
    theta        = theta,
    edf_prof     = edf_prof,
    edf_q50      = edf_q50,
    edf_mo       = edf_mo,
    k_index_prof = k_index_prof,
    kp_prof      = kp_prof,
    k_index_q50  = k_index_q50,
    kp_q50       = kp_q50,
    k_index_mo   = k_index_mo,
    kp_mo        = kp_mo
  )
}


modelos_guardados <- list()
resumen_list      <- list()

for (spp in especies_modelo) {
  
  datos_spp <- lancesspp_modelo %>%
    filter(especie == spp)
  
  for (k in k_grid) {
    
    cat(sprintf("Ajustando %s, k = %d\n", spp, k))
    
    mod  <- ajustar_modelo3_k(datos_spp, k)
    info <- extraer_info_modelo(mod) %>%
      mutate(
        especie   = spp,
        k         = k,
        n_obs     = nrow(datos_spp),
        pct_ceros = mean(datos_spp$densidad_km2 == 0) * 100
      ) %>%
      relocate(especie, k, n_obs, pct_ceros)
    
    modelos_guardados[[paste(spp, k, sep = "_k")]] <- mod
    resumen_list[[paste(spp, k, sep = "_k")]]      <- info
  }
}

resumen_df <- bind_rows(resumen_list) %>%
  arrange(especie, k)

print(resumen_df)

write.csv(
  resumen_df,
  "08_tuning_k_modelo3_resumen.csv",
  row.names = FALSE
)

# Seleccionar mejor k por especie (mínimo AIC)
mejor_k <- resumen_df %>%
  group_by(especie) %>%
  slice_min(order_by = AIC, n = 1, with_ties = FALSE) %>%
  ungroup()

print(mejor_k)


write.csv(
  mejor_k,
  "08_tuning_k_modelo3_mejor_k_porespecie.csv",
  row.names = FALSE
)


# B) Modelos GAM Tweedie con select = TRUE
resultados_gam_select <- list()

for (spp in especies_modelo) {
  cat("\nEspecie:", spp, "(select = TRUE)\n")
  
  datos_spp <- lancesspp_modelo %>%
    filter(especie == spp)
  
  modelo_gam_select <- mgcv::gam(
    densidad_km2 ~
      s(profundidad.x, k = 4) +
      s(Q50phi,        k = 4) +
      s(pctMO,         k = 4) +
      zona + ciclo_estacional + exposicion,
    family = mgcv::tw(),
    data   = datos_spp,
    method = "REML",
    select = TRUE
  )
  
  cat(sprintf(" GAM select converged (AIC: %.2f)\n", AIC(modelo_gam_select)))
  
  suma <- summary(modelo_gam_select)
  if (!is.null(suma$s.table) && nrow(suma$s.table) > 0) {
    cat(" Smooth terms (edf = grados de libertad efectivos):\n")
    print(suma$s.table)
  }
  
  resultados_gam_select[[spp]] <- list(
    modelo = modelo_gam_select,
    AIC    = AIC(modelo_gam_select),
    BIC    = BIC(modelo_gam_select)
  )
  
  cat("\n")
}

# ----------------------------------------------------------------------
# Guardar checkpoints de tuning
# ----------------------------------------------------------------------

if (!dir.exists("checkpoints")) dir.create("checkpoints")

saveRDS(
  list(
    resumen_k  = resumen_df,
    modelos_k  = modelos_guardados,
    gam_select = resultados_gam_select
  ),
  "checkpoints/tuning_k_modelo3.rds"
)

cat("Archivos generados:\n")
cat(" - 08_tuning_k_modelo3_resumen.csv\n")
cat(" - 08_tuning_k_modelo3_mejor_k_porespecie.csv\n")
cat(" - checkpoints/tuning_k_modelo3.rds\n")

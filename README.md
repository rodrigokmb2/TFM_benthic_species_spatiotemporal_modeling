# Modelización espaciotemporal de especies bentónicas en ecosistemas costeros de Galicia

Este repositorio contiene los scripts desarrollados en R utilizados en los análisis
presentados en el Trabajo Fin de Máster:

**“Modelización espaciotemporal de especies bentónicas en ecosistemas costeros de Galicia”**

Autor: Rodrigo Alba Salgueiro  
Programa: Máster Universitario en Bioinformática y Bioestadística  

---

## Descripción general

El presente repositorio acompaña al Trabajo Fin de Máster y reúne los scripts
empleados para el análisis de la distribución espaciotemporal de especies bentónicas
en ecosistemas costeros de Galicia. El estudio integra datos biológicos procedentes de
campañas científicas de arrastre con predictores ambientales de naturaleza física y
sedimentaria, aplicando modelos de distribución de especies basados en enfoques
estadísticos avanzados.

---

## Estructura del repositorio

- `scripts/`: Scripts en R organizados según las distintas fases del flujo de análisis.
- `data/`: Datos de entrada utilizados en el estudio (no disponibles públicamente).
  Véase `data/README_data.txt` para la descripción de su estructura y variables.

---

## Flujo de análisis

Los scripts están diseñados para ejecutarse de forma secuencial, siguiendo el flujo
de trabajo analítico del TFM:

1. `00_preparacion_datos.R`  
2. `01_analisis_comunidad.R`  
3. `02_exploracion_covariables.R`  
4. `03_preparacion_dataset_modelos.R`  
5. `04_seleccion_modelos.R`  
6. `05_tuning_gam_tweedie.R`  

Cada script corresponde a una fase específica del análisis y genera objetos intermedios
utilizados en las etapas posteriores. Una descripción detallada de cada script se
encuentra disponible en `scripts/README_scripts.txt`.

---

## Software

Los análisis se realizaron utilizando R (versión ≥ 4.5.0) y RStudio como entorno de
desarrollo integrado.

---

## Autor

Rodrigo Alba Salgueiro  
GitHub: https://github.com/rodrigokmb2

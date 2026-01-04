# TFM_benthic_species_spatiotemporal_modeling
R scripts for the spatiotemporal modelling of benthic species in coastal ecosystems of Galicia

# Spatiotemporal modelling of benthic species in coastal ecosystems of Galicia

This repository contains the R scripts used in the analyses presented in the Master's Thesis:

**"Modelización espaciotemporal de especies bentónicas en ecosistemas costeros de Galicia"**

## Repository structure

- `scripts/`: R scripts corresponding to each phase of the analysis workflow.
- `data/`: Input data (not publicly available). See `README_data.txt` for details.

## Analysis workflow

Scripts should be executed in the following order:

1. `00_preparacion_datos.R`
2. `01_analisis_comunidad.R`
3. `02_exploracion_covariables.R`
4. `03_preparacion_dataset_modelos.R`
5. `04_seleccion_modelos.R`
6. `05_tuning_gam_tweedie.R`

## Software

Analyses were performed using R (>= 4.5.0) and RStudio.

## Author

Rodrigo Alba Salgueiro

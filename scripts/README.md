Scripts used in the Master's Thesis analysis workflow.

SCRIPTS DESCRIPTION

This directory contains the R scripts used in the analyses presented in the Master's Thesis:

"Modelización espaciotemporal de especies bentónicas en ecosistemas costeros de Galicia"

The scripts are organised following the logical sequence of the analysis workflow.

------------------------------------------------------------
00_preparacion_datos.R
------------------------------------------------------------
Purpose:
- Import raw biological, sedimentary and haul-level data.
- Perform data cleaning, validation and formatting.
- Generate derived variables (e.g. density per km²).
Outputs:
- Cleaned and harmonised datasets used in subsequent analyses.

------------------------------------------------------------
01_analisis_comunidad.R
------------------------------------------------------------
Purpose:
- Characterise benthic community structure.
- Compute alpha diversity indices (Shannon, Simpson, richness, Pielou).
- Perform multivariate analyses (nMDS, SIMPER).
Outputs:
- Community-level summary tables and figures.

------------------------------------------------------------
02_exploracion_covariables.R
------------------------------------------------------------
Purpose:
- Explore relationships among environmental covariates.
- Assess collinearity and variable distributions.
- Select candidate predictors for modelling.
Outputs:
- Correlation matrices and exploratory plots.

------------------------------------------------------------
03_preparacion_dataset_modelos.R
------------------------------------------------------------
Purpose:
- Merge biological, environmental and spatial data.
- Prepare final species-specific modelling datasets.
Outputs:
- Model-ready datasets.

------------------------------------------------------------
04_seleccion_modelos.R
------------------------------------------------------------
Purpose:
- Fit alternative species distribution models (GLM and GAM).
- Compare models using AIC, BIC and explained deviance.
Outputs:
- Model comparison tables and selected final models.

------------------------------------------------------------
05_tuning_gam_tweedie.R
------------------------------------------------------------
Purpose:
- Fine-tune GAM Tweedie models.
- Evaluate sensitivity to smoothing parameters and power estimation.
Outputs:
- Final model diagnostics and effect plots.

------------------------------------------------------------
GENERAL NOTES
------------------------------------------------------------

- Scripts are intended to be executed sequentially.
- Input data are not publicly available; see `data/README_data.txt`.
- All scripts include internal text describing inputs and outputs.

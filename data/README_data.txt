DATA AVAILABILITY AND DESCRIPTION

The data used in this study are not publicly available due to ongoing and planned scientific publications by the research groups and investigators responsible for data acquisition and curation. Access to the raw datasets is therefore restricted at this stage.

This directory describes the structure, content, and role of the datasets used in the analyses presented in the Master's Thesis:

"Modelización espaciotemporal de especies bentónicas en ecosistemas costeros de Galicia"

The datasets were generated within scientific bottom-trawl surveys conducted in the Gulf of Artábro (NW Iberian Peninsula) and are part of institutional research projects coordinated by the Instituto Español de Oceanografía (IEO–CSIC).

----------------------------------------------------------------------
1. datos_lances.csv
----------------------------------------------------------------------

Description:
Metadata associated with each bottom-trawl haul (lance). This table provides spatial, temporal, and operational information used to define the sampling effort and to compute species densities.

Main variables:
- lance: Unique identifier of the trawl haul.
- estacion: Sampling station identifier.
- mes: Month of sampling.
- lat_firmes / lon_firmes: Latitude and longitude at the start of the haul.
- lat_virado / long_virado: Latitude and longitude at the end of the haul.
- area_arrastrada: Estimated swept area (km²).
- zona: Geographic sector (e.g., A Coruña, Ares–Betanzos, Orzán–Riazor).
- exposicion: Binary variable describing sediment exposure (low / high).
- sustrato: Dominant substrate category.
- profundidad: Bottom depth (m).
- estacion.1: Seasonal classification (warm / cold season).

Role in analysis:
Used to define the spatial and temporal structure of the sampling design and to calculate standardized species densities (individuals·km⁻²).

----------------------------------------------------------------------
2. datos_especies.csv
----------------------------------------------------------------------

Description:
Species-level biological data collected for each trawl haul. This table contains abundance information for benthic taxa captured during the surveys.

Main variables:
- estacion: Sampling station identifier.
- lance: Trawl haul identifier.
- cod: Taxonomic code.
- comercial: Indicator of commercial interest.
- nombre_comun: Common name (Spanish).
- nombre_Galicia: Common name (Galician).
- FAO_code: FAO species code.
- especie: Scientific name.
- numero_ejemplares: Number of individuals captured in the haul.

Role in analysis:
Used to characterise community composition, calculate diversity indices, and derive species-specific density metrics for subsequent modelling.

----------------------------------------------------------------------
3. datos_sedimento.csv
----------------------------------------------------------------------

Description:
Sedimentological and environmental characterisation of the sampling stations. This table contains physical and granulometric properties of the seabed.

Main variables:
- estacion: Sampling station identifier.
- zona: Geographic sector.
- profundidad: Bottom depth (m).
- Posicion_Sed: Sediment position within the station.
- pctMO: Percentage of organic matter.
- pctArenaGruesa: Percentage of coarse sand.
- pctArenaFina: Percentage of fine sand.
- pctLodo: Percentage of mud.
- Q50µm: Median grain size (micrometres).
- Q50phi: Median grain size (phi units).
- So: Sorting coefficient.

Role in analysis:
Used as environmental predictors in species distribution models and to explore sediment–community relationships.

----------------------------------------------------------------------
REPRODUCIBILITY NOTE
----------------------------------------------------------------------

Although the original datasets cannot be redistributed, all R scripts provided in the repository are fully documented and structured to allow reproducibility using equivalent datasets with the same variable structure. Variable names and expected formats correspond exactly to those described above.

For access to the original data, interested researchers should contact the corresponding research groups at IEO–CSIC, subject to data availability and institutional policies.

----------------------------------------------------------------------
AUTHOR
----------------------------------------------------------------------

Rodrigo Alba Salgueiro
Master's Thesis – MSc in Bioinformatics and Biostatistics

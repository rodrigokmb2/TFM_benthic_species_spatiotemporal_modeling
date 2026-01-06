DESCRIPCIÓN DE LOS SCRIPTS

Este directorio contiene los scripts desarrollados en R utilizados en los análisis
presentados en el Trabajo Fin de Máster:

"Modelización espaciotemporal de especies bentónicas en ecosistemas costeros de Galicia"

Los scripts están organizados siguiendo la secuencia lógica del flujo de trabajo
analítico, desde la preparación de los datos hasta la selección y ajuste de los
modelos estadísticos.

------------------------------------------------------------
00_preparacion_datos.R
------------------------------------------------------------
Objetivo:
- Importar los datos biológicos, sedimentarios y de lances de arrastre.
- Realizar procesos de depuración, validación y formateo de los datos.
- Generar variables derivadas, como la densidad estandarizada (individuos·km⁻²).
- Salidas: Conjuntos de datos depurados y armonizados utilizados en los análisis posteriores.

------------------------------------------------------------
01_analisis_comunidad.R
------------------------------------------------------------
Objetivo:
- Caracterizar la estructura de la comunidad bentónica.
- Calcular índices de diversidad alfa (Shannon, Simpson, riqueza y equitatividad de Pielou).
- Realizar análisis multivariantes (nMDS, SIMPER).
- Salidas: Tablas resumen y figuras a nivel de comunidad.

------------------------------------------------------------
02_exploracion_covariables.R
------------------------------------------------------------
Objetivo:
- Explorar las relaciones entre las covariables ambientales.
- Evaluar la colinealidad y la distribución de las variables predictoras.
- Seleccionar covariables candidatas para la modelización.
- Salidas: Matrices de correlación y gráficos exploratorios.

------------------------------------------------------------
03_preparacion_dataset_modelos.R
------------------------------------------------------------
Objetivo:
- Integrar la información biológica, ambiental y espacial.
- Preparar los conjuntos de datos finales específicos por especie para la modelización.
- Salidas: Conjuntos de datos listos para la aplicación de modelos de distribución de especies.

------------------------------------------------------------
04_seleccion_modelos.R
------------------------------------------------------------
Objetivo:
- Ajustar modelos alternativos de distribución de especies (GLM y GAM).
- Comparar modelos mediante criterios de información (AIC, BIC) y desviación explicada.
- Salidas: Tablas de comparación de modelos y selección de los modelos finales.

------------------------------------------------------------
05_tuning_gam_tweedie.R
------------------------------------------------------------
Objetivo:
- Ajustar y refinar los modelos GAM con distribución Tweedie.
- Evaluar la sensibilidad de los modelos a los parámetros de suavizado (k) y a la estimación
  del parámetro de potencia (p).
  Salidas: Diagnósticos finales de los modelos y gráficos de efectos parciales.

------------------------------------------------------------
NOTAS GENERALES
------------------------------------------------------------

- Los scripts están diseñados para ejecutarse de forma secuencial.
- Los datos de entrada no son de acceso público; véase `data/README_data.txt` para
  la descripción de su estructura y variables.
- Todos los scripts incluyen documentación interna que describe sus entradas,
  salidas y procedimientos analíticos.

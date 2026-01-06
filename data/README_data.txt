DISPONIBILIDAD Y DESCRIPCIÓN DE LOS DATOS

Los datos utilizados en este estudio no son de acceso público, ya que forman parte de
publicaciones científicas en curso y planificadas por los grupos de investigación e
investigadores responsables de su adquisición y curación. En consecuencia, el acceso
a los conjuntos de datos originales se encuentra restringido en esta fase.

Este directorio describe la estructura, el contenido y el papel de los conjuntos de
datos empleados en los análisis presentados en el Trabajo Fin de Máster:

"Modelización espaciotemporal de especies bentónicas en ecosistemas costeros de Galicia"

Los datos fueron generados en el marco de campañas científicas de arrastre de fondo
realizadas en el Golfo Ártabro (noroeste de la Península Ibérica) y forman parte de
proyectos de investigación institucionales coordinados por el Instituto Español de
Oceanografía (IEO–CSIC), como por ejemplo, el proyecto BIGA (Biodiversidad del Golfo Ártabro).

----------------------------------------------------------------------
1. datos_lances.csv
----------------------------------------------------------------------

Descripción:
Metadatos asociados a cada lance de arrastre de fondo. Esta tabla proporciona
información espacial, temporal y operativa utilizada para definir el esfuerzo de
muestreo y para el cálculo de densidades estandarizadas de las especies.

Variables principales:
- lance: Identificador único del lance de arrastre.
- estacion: Identificador de la estación de muestreo.
- mes: Mes de realización del muestreo.
- lat_firmes / lon_firmes: Latitud y longitud al inicio del lance.
- lat_virado / long_virado: Latitud y longitud al final del lance.
- area_arrastrada: Área barrida estimada (km²).
- zona: Sector geográfico (A Coruña, Ares–Betanzos, Orzán–Riazor).
- exposicion: Variable binaria que describe la perturbación física combinada del sedimento del fondo marino, resultante de la acción conjunta de mareas, oleaje y corrientes inducidas por el viento (0 = baja / 1 = alta).
- sustrato: Categoría dominante del sustrato (3 categorías: arena fina, arena gruesa, grava o cascajos).
- ciclo_estacional: Clasificación del ciclo estacional (estación cálida, de abril a septiembre / estación fría, de octubre a marzo).
- profundidad: Profundidad del fondo (m).

Papel en el análisis:
Utilizado para definir la estructura espacial y temporal del diseño de muestreo y para
el cálculo de densidades estandarizadas (individuos·km⁻²).

----------------------------------------------------------------------
2. datos_especies.csv
----------------------------------------------------------------------

Descripción:
Datos biológicos a nivel de especie recopilados para cada lance de arrastre. Esta tabla
contiene información sobre la abundancia de los taxones bentónicos capturados durante
las campañas.

Variables principales:
- estacion: Identificador de la estación de muestreo.
- lance: Identificador del lance de arrastre.
- cod: Código taxonómico.
- comercial: Indicador de interés comercial.
- nombre_comun: Nombre común en castellano.
- nombre_Galicia: Nombre común en gallego.
- FAO_code: Código FAO de la especie.
- especie: Nombre científico.
- numero_ejemplares: Número de individuos capturados en el lance.

Papel en el análisis:
Utilizado para caracterizar la composición de la comunidad bentónica, calcular índices
de diversidad y derivar métricas de densidad específicas por especie para su posterior
modelización.

----------------------------------------------------------------------
3. datos_sedimento.csv
----------------------------------------------------------------------

Descripción:
Caracterización sedimentológica y ambiental de las estaciones de muestreo. Esta tabla
incluye propiedades físicas y granulométricas del fondo marino.

Variables principales:
- estacion: Identificador de la estación de muestreo.
- zona: Sector geográfico (A Coruña, Ares–Betanzos, Orzán–Riazor).
- profundidad: Profundidad del fondo (m).
- Posicion_Sed: Posición de la muestra de sedimento dentro de la estación.
- pctMO: Porcentaje de materia orgánica (contenido orgánico total respecto al peso seco de la muestra).
- pctArenaGruesa: Porcentaje de arena gruesa.
- pctArenaFina: Porcentaje de arena fina.
- pctLodo: Porcentaje de lodo (suma de limos + arcillas, fracción más fina del sedimento).
- Q50µm: Tamaño medio del grano (micrómetros). Mediana del tamaño de grano expresada en micrómetros
- Q50phi: Tamaño medio del grano (unidades phi). La misma mediana de tamaño de grano pero transformada a unidades phi 
- So: Coeficiente de selección (sorting). Cuantifica cuán homogénea es la distribución de tamaños de grano

Papel en el análisis:
Utilizado como conjunto de predictores ambientales en los modelos de distribución de
especies y para explorar las relaciones entre el sedimento y la estructura de la
comunidad bentónica.

----------------------------------------------------------------------
NOTA SOBRE REPRODUCIBILIDAD
----------------------------------------------------------------------

Aunque los conjuntos de datos originales no pueden ser redistribuidos, todos los
scripts en R proporcionados en este repositorio están documentados y estructurados de
forma que permiten reproducir el flujo completo de análisis utilizando datos
equivalentes con la misma estructura y nomenclatura de variables. Los nombres de las
variables y los formatos esperados se corresponden exactamente con los descritos en
este documento.

Para solicitar acceso a los datos originales, los investigadores interesados deberán
contactar con los grupos de investigación correspondientes del IEO–CSIC, de acuerdo
con las políticas institucionales y la disponibilidad de los datos.

----------------------------------------------------------------------
AUTOR
----------------------------------------------------------------------

Rodrigo Alba Salgueiro  
Trabajo Fin de Máster – Máster Universitario en Bioinformática y Bioestadística

# COVID_Baq
Estimación de curvas epidémicas de COVID 19 en Barranquilla, localidades y otras ciudades en Colombia. 

Ultima actualización: 
Noviembre 19 de 2020

Responsable del código y a quien deben ser dirigidas las inquietudes:
Juan Pablo Gómez
echeverrip@uninorte.edu.co

El código contenido en esta carpeta es el requerido para recrear la sección de predicciones y calcular el Rt en el aplicativo web entregado como parte del documento base que acompaña el presente informe.

Cada archivo contiene diferentes elementos del código fuente el cual ha sido debidamente anotado. En este archivo se detallara el contenido de cada uno de los códigos fuente.

El paquete depende de las siguientes librerías y sus dependencias:

1. EnvStats
2. mvtnorm
3. rjags
4. parallel
5. dclone
6. MASS
7. deSolve

Todas las librerías son llamadas desde el archivo COVD_Baq_librerias.R asi como los archivos del código fuente. Por consiguiente este archivo será el primero que debe ser leido como código fuente. Al cargarlo en R, le preguntara la ruta para buscar el resto del código fuente.

La estimación esta basada en cinco archivos principales.

1. COVID_Baq_organizador_datos.R
Funciones para descargar y organizar los datos de la página oficial del INS. Adicionalmente hay una opción para organizar los datos de la Alcaldía el cual contiene información sobre la localidad de cada uno de los casos. En la función se ha anotado el lugar en el que se debe insertar la ruta para el correcto funcionamiento de la función. Ver anotaciones en el código. Adicionalmente, este código asume que la información de la localidad se encuentra guardada en la base de datos en una columna a parte denominada como LOCALIDA y que el resto de columnas y entradas en las columnas son idénticas a la base de datos del INS.

2. COVID_Baq_modelos_logisticos.R
Funciones que contienen los diferentes modelos de crecimiento logistico deterministico y estocástico descritos en el documento base.

3. COVID_Baq_modelos_SIR.R
Funciones con el sistema de ecuaciones diferenciales para un modelo de tipo SIR no estructurado.

4. COVID_Baq_funciones_estimacion.R
Funciones requeridas para estimar los parametros de cada uno de los modelos en COVID_Baq_modelos_logisticos.R

5. COVID_Baq_funciones_prediccion.R
Funciones requeridas para realizar proyecciones a futuro de los modelos logísticos descritos en el documento base.

A partir de estos cinco archivos se construyeron los diferentes resultados obtenidos.

En el archivo COVID_Baq_Estimacion.R contiene el código necesario para realizar la estimación de los modelos generalizados logísticos, las proyecciones para Barranquilla, sus localidades y las demás ciudades. 

En el archivo COVID_Baq_Estimacion_SIR.R se puede encontrar el código necesario para reproducir el modelamiento bajo los modelos compartimentalizados de tipo SIR descritos en el informe adjunto. Tambien es posible recrear la estimación del número de infectados no documentados utilizando el método de Mellis y Littera 2020. 

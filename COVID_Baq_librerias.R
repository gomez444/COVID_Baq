#####################################################################################################
#
# Librerias y código fuente requeridos para correr el paquete de estimación de casos de COVID
# Ultima actualización: 9 de Diciembre de 2020
#
####################################################################################################
# Librerias
library(EnvStats)
library(mvtnorm)
#library(rjags)
library(parallel)
#library(dclone)
library(MASS)
library(deSolve)

# Cargar codigo base requerido 
#path	<-	readline(prompt="Defina la ruta para buscar los archivos con las funciones: \n")
source(paste(path,'COVID_Baq_organizador_datos.R',sep=''), chdir = TRUE,encoding = 'UTF-8')
source(paste(path,'COVID_Baq_modelos_logisticos.R',sep=''), chdir = TRUE,encoding = 'UTF-8')
source(paste(path,'COVID_Baq_funciones_estimacion.R',sep=''), chdir = TRUE,encoding = 'UTF-8')
source(paste(path,'COVID_Baq_funciones_prediccion.R',sep=''), chdir = TRUE,encoding = 'UTF-8')
source(paste(path,'COVID_Baq_modelos_SIR.R',sep=''), chdir = TRUE,encoding = 'UTF-8')

# Función para obtener la paleta de colores
my.col.ramp	<-	colorRampPalette(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','grey','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))
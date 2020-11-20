#####################################################################################################
#
# Funciones para evaluar y calcular la verosimilitud de diferentes modelos de crecimiento logístico 
# Ultima actualización: 18 de Noviembre de 2020
#
####################################################################################################
##################################################
#
# Modelos deterministicos
#
##################################################
# modelo Gompertz deterministico
# b0 =  valor del intercepto para el cambio de la tasa de crecimiento en el tiempo
# b1 = valor de la pendiente para el cambio de la tasa de crecimiento en el tiempo
# c = parámetro que controla la infección dependiente de la densidad
# no = valor inicial para la estimación
# len = longitud total del número de días para estimar el modelo
# to = valor inicial del tiempo en dias desde el inicio de la epidemia
# Salida
# Vector con el número casos de longitud len
gompertz.det	<- function(b0,b1,c,no,len,to){
		
	Nt<-rep(0,len)
	Nt[1]<-no
	for(i in 2:len){
		a	<-	exp(b0+b1*(to+(i-1)))
		Nt[i]	<-	exp(a+c*log(Nt[(i-1)]))
		
	}
	return(Nt)
}

####################################################################################################
# Modelo Gompertz determinístico condicional en las observaciones
# b0 =  valor del intercepto para el cambio de la tasa de crecimiento en el tiempo
# b1 = valor de la pendiente para el cambio de la tasa de crecimiento en el tiempo
# c = parámetro que controla la infección dependiente de la densidad
# no = valor inicial para la estimación
# len = longitud total del número de días para estimar el modelo
# to = valor inicial del tiempo en dias desde el inicio de la epidemia
# dat = vector de observaciones de número acumulado de infectados a traves del tiempo
# Salida
# Vector con el número de casos de longitud len

gompertz.cond	<-	function(b0,b1,c,len,to,dat){
	Nt<-rep(0,len)
	Nt[1]<-dat[1]
	for(i in 2:len){
		a	<-	exp(b0+b1*(to+(i-1)))
		Nt[i]	<-	exp(a+c*log(dat[(i-1)]))
	}
	return(Nt)
}

####################################################################################################
# Modelo logístico Generalizado (Richards function)
# K = Capacidad de carga o número total de infectados al final de la epidemia.
# r = tasa de crecimiento
# theta = tiempo desde el inicio de la epidemia en el que se alcanza el punto de inflexión
# alpha = parámetro que controla la asimetria en la curva de crecimiento
# t = vector de tiempo en dias en los que será evaluado el modelo. 
# Salida
# Vector con el número de casos de longitud t

gen.log.mod	<-	function(K,r,theta,alpha,t){
	
	K*(1+alpha*exp(-r*(t-theta)))^(-1/alpha)
	
}

####################################################################################################
# Modelo de crecimiento propuesto por con Mellis y Litera.
# c1 = parametro 1.
# c2 = parametro 2.
# c3 = parametro 3.
# x = vector de tiempo para calcular el número predicho de eventos.
recov.t	<-	function(c1,c2,c3,x){
	
	c1*(tanh(c2*x-c3)+tanh(c3))
	
}


#################################################################################### 
# Modelos para estimación del modelo Gompertz estocastico utilizando
# clonación de datos
# Modelo estocástico de Gompertz con tasa de crecimiento constante (a)
# Modelo en lenguaje jags para se evaluado por rjags

const.growth.model.dc <- function(){
  
	lsig ~ dnorm(0,0.0001)
  	sig <- exp(lsig)
	sigsq.inv <- pow(sig,-2)
  	
  	cc ~	dunif(0,1)
    loga ~ dnorm(0, 0.001)
    a    <- exp(loga)

    for(k in 1:K){
      Xt[1,k]	 ~ dnorm(a+(cc*log(nos)),sigsq.inv)
      Nt[1,k] <- exp(Xt[1,k])
      N.obs[1,k] ~ dpois(Nt[1,k])#;T(nos,)
      for(j in 2:lens){ 
          Xt[j,k] ~ dnorm((a + (cc)*Xt[(j-1),k]), sigsq.inv)
          Nt[j,k]	<-	exp(Xt[j,k])
          N.obs[j,k] ~ dpois(Nt[j,k])#;T(Nt[j-1,k],)
      }
    }
  }  

####################################################################################
# Modelo estocástico de Gompertz con tasa de crecimiento variable (a)
# Modelo en lenguaje jags para se evaluado por rjags

var.growth.model.dc <- function(){

  beta0	~	dnorm(0,0.001)
  beta1 ~ 	dnorm(0,0.001)
     
  lsig ~ dnorm(0,0.0001)
  sig <- exp(lsig)
  sigsq.inv <- pow(sig,-2)
  
  cc ~ dunif(0,1)

   for(k in 1:K){
		a[1,k]	<-	exp(beta0+beta1*1)
    		Xt[1,k]	 ~ dnorm(a[1,k]+(cc*log(nos)),sigsq.inv)
    		Nt[1,k] <- exp(Xt[1,k])
    		
      for(j in 2:lens){
      	  a[j,k]	<-	exp(beta0+beta1*j)
          Xt[j,k] ~ dnorm((a[j,k] + cc*Xt[(j-1),k]), sigsq.inv)
          Nt[j,k]	<-	exp(Xt[j,k])
      }
    }
    for(k in 1:K){
			N.obs[1,k]~dpois(Nt[1,k])#;T(nos,)
    		for(j in 2:lens){
    			N.obs[j,k]~dpois(Nt[j,k])#;T(Nt[j-1,k],)
    		}
    }
}  

####################################################################################
#####  Estimadores de Kalman para los modelos Gompertz estocásticos
# Modelo estocástico de Gompertz con tasa de crecimiento constante (a)
# Modelo en lenguaje jags para se evaluado por rjags

const.growth.kalman <- function(){
  
    Xt[1]	 ~ dnorm(a+(cc*log(nos)),sigsq.inv)
    Nt[1] <- exp(Xt[1])
	N.obs[1] ~ dpois(Nt[1])
    for(j in 2:lens){
          Xt[j] 	~ dnorm((a + cc*Xt[(j-1)]), sigsq.inv)
          Nt[j]	<-	exp(Xt[j])
		  N.obs[j] ~ dpois(Nt[j])
    }
}  

####################################################################################
# Modelo estocástico de Gompertz con tasa de crecimiento variable (a)
# Modelo en lenguaje jags para se evaluado por rjags

var.growth.kalman <- function(){
  	
  	a[1]	<-	exp(beta0+beta1*1)
    Xt[1]	 ~ dnorm(a[1]+(cc*log(nos)),sigsq.inv)
    Nt[1] <- exp(Xt[1])
    N.obs[1] ~ dpois(Nt[1])#;T(nos,)
    
    for(j in 2:lens){
    		a[j]	 	<- exp(beta0+beta1*j)	
        Xt[j] ~ dnorm((a[j] + (cc)*Xt[(j-1)]), sigsq.inv)
        Nt[j]	<-	exp(Xt[j])
        N.obs[j]~ dpois(Nt[j])#;T(Nt[j-1],)
    }
}  

####################################################################################
# Funciones para la simulación de trayectorias bajo los modelos de Gompertz estocásticos
# Modelo estocástico de Gompertz con tasa de crecimiento constante (a)
# parms = vector de longitud 3 con los estimadores de máxima verosimilitud del modelo
# gompertz estocastico de crecimiento constante.
# len = número de pasos para la simulación
# no = valor incial de la simulación
# Salida:
# Vector de longitud len con el resultado de la simulación dados los parámetros

const.growth.sim <- function(parms, len, no){
  
  a      <- parms[1]
  cc     <- parms[2]
  sig    <- parms[3]
  
  Nt    <- rep(0,len)
  Nt[1] <- no
  
  for(i in 2:len){
    Xt   	<- rnorm(n=1,a+(cc)*log(Nt[(i-1)]), sd=sig)
    Nt[i]	<- rpois(1,exp(Xt))
  }
  return(Nt)
}

####################################################################################
# Modelo estocástico de Gompertz con tasa de crecimiento variable (a)
# parms = vector de longitud 3 con los estimadores de máxima verosimilitud del modelo
# gompertz estocastico de crecimiento variable.
# parms = vector de longitud 3 con los estimadores de máxima verosimilitud del modelo
# gompertz estocastico de crecimiento constante.
# len = número de pasos para la simulación
# no = valor incial de la simulación
# Salida:
# Vector de longitud len con el resultado de la simulación dados los parámetros

var.growth.sim <- function(parms, len, no,to){
  
  cc	     <- parms[1]
  beta0	 <- parms[2]	
  beta1	 <-	parms[3]
  sig    <- parms[4]
  
  Nt    <- rep(0,len)
  Nt[1] <- no
  
  for(i in 2:len){
	a		<-	exp(beta0+beta1*(to+i-1))
    Xt   	<- rnorm(n=1,a+(cc)*log(Nt[(i-1)]), sd=sig)
	Nt[i]	<-	max(Nt[i-1],rpois(1,exp(Xt)))
  
  }
  return(Nt)
}

####################################################################################
# Modelo de Gompertz deterministico condicional
# b0 =  valor del intercepto para el cambio de la tasa de crecimiento en el tiempo
# b1 = valor de la pendiente para el cambio de la tasa de crecimiento en el tiempo
# c = parámetro que controla la infección dependiente de la densidad
# no = valor inicial para la estimación
# len = longitud total del número de días para estimar el modelo
# to = valor inicial del tiempo en dias desde el inicio de la epidemia
# dat = vector de observaciones de número acumulado de infectados a traves del tiempo
# Salida:
# Vector de longitud len con el resultado de la simulación dados los parámetros

gompertz.cond.sim	<-	function(b0,b1,cc,len,to,dat){
	Nt<-rep(0,len)
	Nt[1]<-dat[1]
	for(i in 2:len){
		a	<-	exp(b0+b1*(to+(i-1)))
		Nt[i]	<-	max(Nt[i-1],rpois(1,exp(a+cc*log(dat[(i-1)]))))
	}
	return(Nt)
}
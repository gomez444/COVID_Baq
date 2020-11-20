#### Modelos SIR
# Los siguientes modelos son sistemas de ecuaciones diferenciales 
# de la forma requerido por el paquete deSolve para encontrar la solución del sistema.
# Sistema de ecuaciones diferenciales SIR simple sin estructura.
# t =  vector con la serie de tiempo para el cual se quiere evaluar el sistema.
# state =  vector con valores iniciales para los estados S,I,R,C,D
# parameters = valores de b0, beta, theta, gamma y mu para los cuales se quiere la solución del sistema
# b0  = Coeficiente que determina la velocidad del cambio en la proporción de susceptibles disponibles para contagio
# beta =  tasa de contagio.
# theta = Número de individuos immigrantes contagiados.
# gamma = tasa de recuperación de infectados
# mu = tasa de mortalidad de infectados.
# Simple Model
sir.simple <-function(t, state, parameters) {
	with(as.list(c(state, parameters)),{
		quar.prop	<-	exp(-b0*t)
		theta<-theta.func(t,theta)
		dS 	<- -(beta*(I/(S+I+R))*S) - quar.prop*S
		dI 	<- (beta*(I/((S+I+R)))*S) - (gamma + mu)*I + theta
		dR	<-	(gamma*I)
		dC	<- (beta*(I/((S+I+R)))*S)+theta
		dD 	<-  mu*I
		# return the rate of change
		list(c(dS,dI,dR,dC,dD))
	}) # end with(as.list ...
}
##########################################################################################
##########################################################################################
##########################################################################################
##########
##########
########## Modelamiento de COVID con base en modelos SIR: Ver descripción del modelo y 
########## supuestos en el documento adjunto.
##########
##########################################################################################
##########################################################################################
##########################################################################################

####################################################################################################
#
# Ultima actualización: 18 de Noviembre de 2020
#
####################################################################################################

# Obtener los datos
path.dat		<-	" " #definir el directorio para guardar o buscar los datos. 
download		<-	FALSE #Definir si se quiere descargar los datos o buscar en el directorio.
path.figs		<-	" "# Ruta para guardar las figuras

data.ciudad	<-	getCOVID.data(wd=path.dat,location="BARRANQUILLA",type="Municipio"
								,download=download,date.type="Diagnostico",data.origin="INS")

# Estimación de parametros
# Datos fijos del modelo
pop.file		<-	read.csv(paste(path.dat,"edades.csv",sep=""))
pop	<-	sum(pop.file$value[pop.file$Localidad=="consolidado"])
air.close.date	<-	as.Date("2020-03-23") # Fecha de cierre de fronteras
first.date	<-	as.Date("2020-03-16") # Fecha de reporte de primer caso
air.close	<-	as.numeric(air.close.date-first.date)
dat	<-	t(data.ciudad$BARRANQUILLA[data.ciudad$BARRANQUILLA[,1]>0,c(1,3,5)])
new.cases	<-	data.ciudad$BARRANQUILLA[data.ciudad$BARRANQUILLA[,1]>0,2]
guess	<-	numeric(3) # Número de parametros para estimar dentro del modelo SIR
gamma	<-	1/14 # Tasa de recuperación. En promedio dos semanas de infecciosidad
activos	<-	dat[1,]-dat[3,]
muertos.dia	<-	data.ciudad$Barranquilla[data.ciudad$Barranquilla[,1]>0,4]
mort.diaria<-muertos.dia[17:length(activos)]/activos[1:(length(activos)-16)]
q95		<-	quantile(mort.diaria,prob=0.95)
q5		<-	quantile(mort.diaria,prob=0.05)	
mu.obs	<-	mean(mort.diaria[mort.diaria<q95])
mu		<-	0.0034# # Tasa de mortalidad.

# Crear una lista con los datos para la optimización.
data4optim		<-	list(ts=dat,mu=mu,gamma=gamma,pop=pop)

# Funciones del modelo. 
# Permitimos inmigración de individuos infectados a la población 
# hasta el momento en el que se cierran fronteras
theta.func	<-	function(t,theta){if(t<=air.close){theta}else{0}}

# Optimización de los parametros beta (tasa de contagio), beta0 y beta1 (coeficientes de cambio en la población)
# de susceptibles 
ith.optim		<-	optim(par=guess,fn=sir.simp.negloglik,data=data4optim
							,method="Nelder-Mead",hessian=TRUE)
obs.len	<-	ncol(dat)
pars.se	<-	sqrt(diag(ginv(ith.optim$hessian)))
pars.ci	<-	matrix(c(exp(ith.optim$par[1:2]-1.96*pars.se[1:2]),ith.optim$par[3]-1.96*pars.se[3]
					,exp(ith.optim$par[1:2]+1.96*pars.se[1:2]),ith.optim$par[3]+1.96*pars.se[3])
					,ncol=3,nrow=2,byrow=TRUE)

# Simulación para un año a partir del primer caso.
pars4sim	<-	c(beta=exp(ith.optim$par[1]),mu=mu,gamma=gamma
					,theta=exp(ith.optim$par[2]),b0=ith.optim$par[3])
fechas	<-	seq.Date(first.date,by="days",length.out=obs.len)
fechas.proy<-	seq.Date(first.date,by="days",length.out=365)
sims.state	<-	c(S=pop,I=as.numeric(pars4sim["theta"]),R=0,C=as.numeric(pars4sim["theta"]),D=0)
sims1		<-	ode(y=sims.state,times=1:365,parms=pars4sim,func=sir.simple)

# R0 estimado a partir de la relación entre beta y gamma
R0	<-	pars4sim["beta"]/gamma

# Generar la figura de la estimación de las diferentes series de tiempo
tiff(paste(path.figs,"Ajuste_SIR.tiff",sep=""),width=5,height=6,units="in",res=300
	,type="cairo",compression="lzw",family="times")
par(mfrow=c(3,2),mar=c(3,3,1,1),oma=c(0,0,2,0),tcl=-0.2,mgp=c(1.5,0.5,0))
plot(fechas.proy,sims1[,"C"],ylim=c(0,max(sims1[,"C"])),xlab="Fecha",ylab="Número acumulado de Casos",type="n"
	,cex=0.5,cex.axis=0.7,cex.lab=0.7,bty="l")
points(fechas,dat[1,],pch=16,cex=0.5)
points(fechas.proy,sims1[,"C"],type="l",col="red")

new.cases.pred	<-	c(sims1[1,"C"],sims1[2:365,"C"]-sims1[1:(364),"C"])
plot(fechas.proy,new.cases.pred,type="n",ylab="Número diario de casos nuevos",xlab="Fechas",ylim=c(0,max(new.cases))
	,cex=0.5,cex.axis=0.7,cex.lab=0.7,bty="l")
points(fechas,new.cases,pch=16,cex=0.5)
points(fechas.proy,new.cases.pred,type="l",col="red")

plot(fechas.proy,sims1[,"D"],pch=19,xlab="Fecha",ylab="Número acumulado de muertes",type="n"
	,cex=0.5,cex.axis=0.7,cex.lab=0.7,bty="l")
points(fechas,dat[2,],pch=16,cex=0.5)
points(fechas.proy,sims1[,"D"],type="l",col="red")

plot(fechas.proy,sims1[,"R"],type="n",ylab="Número acumulado de Recuperados",xlab="Fecha"
	,cex=0.5,cex.axis=0.7,cex.lab=0.7,bty="l")
points(fechas,dat[3,],pch=16,cex=0.5)
points(fechas.proy,sims1[,"R"],type="l",col="red")

activos	<-	dat[1,]-dat[2,]-dat[3,]
plot(fechas.proy,sims1[,"I"],type="l",col="red",ylab="Número de Infectados",xlab="Fecha",ylim=c(0,max(activos))
	,cex=0.5,cex.axis=0.7,cex.lab=0.7,bty="l")
points(fechas,activos,pch=16,cex=0.5)

plot(fechas.proy,exp(-pars4sim["b0"]*1:365),col="red",ylab="Fracción de población no infectada removida",type="l",xlab="Fecha"
	,cex=0.5,cex.axis=0.7,cex.lab=0.7,bty="l")
par(fig=c(0,1,0,1),new=TRUE,mar=c(0,3,0,1),oma=c(0,0,0,0))
plot(0,0,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
legend("top",legend=c("Observados","Predichos"),pch=c(16,-1),lty=c(-1,1),col=c("black","red"),bty="n",horiz=TRUE)
dev.off()

###############################################################################################
#
#
# Simulación de diferentes escenarios para un año a partir del incio de la pandemia
#
#
###############################################################################################
########################################################################
# Escenario 0: Optimista
# La trayectoria continua como termino a fecha del 8 de Septiembre
state.0	<-	sims1[obs.len,-1]
sims.0	<-	ode(y=state.0,times=obs.len:365,parms=pars4sim,func=sir.simple)
sims.0	<-	rbind(sims1,sims.0[-1,])

########################################################################
# Escenario 1:Pesimista 100%
# Toda la población susceptible removida y suceptible al final de la cuarentena que no ha sido infectada entra 
# en la población suceptible 15 días despues de Septiembre 1 cuando hay reactivación total.
sims1		<-	ode(y=sims.state,times=1:(ncol(dat)+7),parms=pars4sim,func=sir.simple)
ntot			<-	pop-tail(sims1[,"I"],n=1)-tail(sims1[,"D"],n=1)-tail(sims1[,"R"],n=1)
state.1.2		<-	sims1[ncol(dat),-1]
state.1.2["S"]	<-	ntot
sims.1.2		<-	ode(y=state.1.2,times=(obs.len+7):365,parms=pars4sim,func=sir.simple)
sims.1			<-	rbind(sims1,sims.1.2[-1,])

########################################################################
# Escenario 1:Pesimista 80%
# Toda la población susceptible removida y suceptible al final de la cuarentena que no ha sido infectada entra 
# en la población suceptible 15 días despues de Septiembre 1 cuando hay reactivación total.
ntot			<-	(pop-tail(sims1[,"I"],n=1)-tail(sims1[,"D"],n=1)-tail(sims1[,"R"],n=1))*0.8
state.1.2["S"]	<-	ntot
sims.1.3		<-	ode(y=state.1.2,times=(obs.len+7):365,parms=pars4sim,func=sir.simple)
sims.1.3		<-	rbind(sims1,sims.1.3[-1,])

########################################################################
# Escenario 1:Pesimista 50%
# Toda la población susceptible removida y suceptible al final de la cuarentena que no ha sido infectada entra 
# en la población suceptible 15 días despues de Septiembre 1 cuando hay reactivación total.
ntot			<-	(pop-tail(sims1[,"I"],n=1)-tail(sims1[,"D"],n=1)-tail(sims1[,"R"],n=1))*0.5
state.1.2["S"]	<-	ntot
sims.1.4		<-	ode(y=state.1.2,times=(obs.len+7):365,parms=pars4sim,func=sir.simple)
sims.1.4		<-	rbind(sims1,sims.1.4[-1,])

########################################################################
# Escenario 1:Pesimista 30%
# Toda la población susceptible removida y suceptible al final de la cuarentena que no ha sido infectada entra 
# en la población suceptible 15 días despues de Septiembre 1 cuando hay reactivación total.
ntot			<-	(pop-tail(sims1[,"I"],n=1)-tail(sims1[,"D"],n=1)-tail(sims1[,"R"],n=1))*0.3
state.1.2["S"]	<-	ntot
sims.1.5		<-	ode(y=state.1.2,times=(obs.len+7):365,parms=pars4sim,func=sir.simple)
sims.1.5		<-	rbind(sims1,sims.1.5[-1,])

########################################################################
# Escenario 1:Pesimista 10%
# Toda la población susceptible removida y suceptible al final de la cuarentena que no ha sido infectada entra 
# en la población suceptible 15 días despues de Septiembre 1 cuando hay reactivación total.
ntot			<-	(pop-tail(sims1[,"I"],n=1)-tail(sims1[,"D"],n=1)-tail(sims1[,"R"],n=1))*0.1
state.1.2["S"]	<-	ntot
sims.1.6		<-	ode(y=state.1.2,times=(obs.len+7):365,parms=pars4sim,func=sir.simple)
sims.1.6		<-	rbind(sims1,sims.1.6[-1,])

fechas	<-	seq.Date(as.Date("2020-03-16"),by="days",length.out=365)

my.cols	<-	my.col.ramp(5)

# Generar la figura para el escenario 1
tiff(paste(path.figs,"Escenario1_SIR.tiff",sep=""),width=5,height=5,units="in",res=300
	,type="cairo",compression="lzw",family="times")
par(mar=c(3,3,1,1),oma=c(0,0,2,0),tcl=-0.2,mgp=c(1.5,0.5,0))
plot(fechas,sims.1[,"I"],type="l",col=my.cols[1],xlab="Fecha",ylab="Número de Infectados")
points(fechas,sims.1.3[,"I"],type="l",col=my.cols[2])
points(fechas,sims.1.4[,"I"],type="l",col=my.cols[3])
points(fechas,sims.1.5[,"I"],type="l",col=my.cols[4])
points(fechas,sims.1.6[,"I"],type="l",col=my.cols[5])
par(fig=c(0,1,0,1),new=TRUE,mar=c(0,3,1,1),oma=c(0,0,0,0))
plot(0,0,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
legend("top",legend=c("100%","80%","50%","30%","10%"),lty=1,col=my.cols,bty="n",horiz=TRUE,cex=0.8)
dev.off()
########################################################################
# Escenario 2: Semi optimista
# La reactivación incrementa el número de individuos susceptibles incrementando en 15% la población suceptible
# cada 2 meses hasta que se alcanza el 100% de la población
sims.2	<-	sims1
for(i in 1:3){
	ith.state	<-	sims.2[nrow(sims.2),-1]
	ntot			<-	(pop-tail(sims.2[,"I"],n=1)-tail(sims.2[,"D"],n=1)-tail(sims.2[,"R"],n=1))*0.15
	ith.state["S"]	<- ntot
	ith.time		<-	tail(sims.2[,"time"],n=1)
	ith.sims		<-	ode(y=ith.state,times=ith.time:(ith.time+60),parms=pars4sim,func=sir.simple)
	sims.2	<-	rbind(sims.2,ith.sims[-1,])
}

sims.2.1	<-	sims1
for(i in 1:3){
	ith.state	<-	sims.2.1[nrow(sims.2.1),-1]
	ntot			<-	(pop-tail(sims.2.1[,"I"],n=1)-tail(sims.2.1[,"D"],n=1)-tail(sims.2.1[,"R"],n=1))*0.10
	ith.state["S"]	<- ntot
	ith.time		<-	tail(sims.2.1[,"time"],n=1)
	ith.sims		<-	ode(y=ith.state,times=ith.time:(ith.time+60),parms=pars4sim,func=sir.simple)
	sims.2.1	<-	rbind(sims.2.1,ith.sims[-1,])
}

sims.2.2	<-	sims1
for(i in 1:3){
	ith.state	<-	sims.2.2[nrow(sims.2.2),-1]
	ntot			<-	(pop-tail(sims.2.2[,"I"],n=1)-tail(sims.2.2[,"D"],n=1)-tail(sims.2.2[,"R"],n=1))*0.05
	ith.state["S"]	<- ntot
	ith.time		<-	tail(sims.2.2[,"time"],n=1)
	ith.sims		<-	ode(y=ith.state,times=ith.time:(ith.time+60),parms=pars4sim,func=sir.simple)
	sims.2.2<-	rbind(sims.2.2,ith.sims[-1,])
}

sims.2.3	<-	sims1
for(i in 1:94){
	ith.state	<-	sims.2.3[nrow(sims.2.3),-1]
	ntot			<-	(pop-tail(sims.2.3[,"I"],n=1)-tail(sims.2.3[,"D"],n=1)-tail(sims.2.3[,"R"],n=1))*0.05
	ith.state["S"]	<- ntot
	ith.time		<-	tail(sims.2.3[,"time"],n=1)
	ith.sims		<-	ode(y=ith.state,times=ith.time:(ith.time+2),parms=pars4sim,func=sir.simple)
	sims.2.3<-	rbind(sims.2.3,ith.sims[-1,])
}

last.date	<-	ncol(dat)
my.cols	<-	my.col.ramp(4)
tiff(paste(path.figs,"Escenario2_SIR.tiff",sep=""),width=5,height=5,units="in",res=300
	,type="cairo",compression="lzw",family="times")
par(mar=c(3,3,1,1),oma=c(0,0,2,0),tcl=-0.2,mgp=c(1.5,0.5,0))
plot(fechas,sims.2[1:365,"I"],type="l",col=my.cols[1],xlab="Fecha",ylab="Número de Infectados")
points(fechas,sims.2.1[1:365,"I"],type="l",col=my.cols[2])
points(fechas,sims.2.2[1:365,"I"],type="l",col=my.cols[3])
points(fechas,sims.2.3[1:365,"I"],type="l",col=my.cols[4])
par(fig=c(0,1,0,1),new=TRUE,mar=c(0,3,0.5,1),oma=c(0,0,0,0))
plot(0,0,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
legend("top",legend=c("15%","10%","5%","5% Constante"),lty=1,col=my.cols,bty="n",ncol=2,cex=0.8)
dev.off()


########################################################################
# Escenario 3: Semi optimista 2
# Incrementar la población a la velocidad que se disminuyo
alpha.t<-function(t,b0,b1){exp(-(b0+b1*(177/t)))}
state.3	<-	c(sims1[obs.len,-1],N=as.numeric(sims1[obs.len,"S"]+sims1[obs.len,"I"]+sims1[obs.len,"R"]))
pars4sim["b1"]	<-	0.05
ntot			<-	pop-tail(sims1[,"I"],n=1)-tail(sims1[,"D"],n=1)-tail(sims1[,"R"],n=1)
pars4sim3	<-	c(pars4sim,Ntot=ntot)
sims.3	<-	ode(state.3,times=177:365,parms=pars4sim3,func=sir.simple.sim)
sims.3	<-	rbind(sims1,sims.3[-1,])

sums	<-	apply(sims.3[,-1],1,sum)
plot(sums)

# Aislamiento ciclico de 15 días por localidades
pop.locs	<-	tapply(pop.file$value,pop.file$Localidad,sum)[-1]

plot(new.cases,,pch=16,ylim=c(0,10000),xlim=c(0,365))
points(sims.1[2:365,"C"]-sims.1[1:364,"C"],type="l",col="blue")
points(sims.0[2:365,"C"]-sims.0[1:364,"C"],type="l",col="red")
points(sims.2[2:364,"C"]-sims.2[1:363,"C"],type="l",col="forestgreen")
points(sims.3[2:364,"C"]-sims.3[1:363,"C"],type="l",col="darkorange")

########################################################################
# Estimación del número de recuperados y el número de infectados no reportados según
# Melis y Littera 2020.
recov.hat	<-	optim(par=rnorm(3),recov.nll,x=dat[3,]+dat[2,],hessian=TRUE)
recov.se	<-	sqrt(ginv(diag(recov.hat$hessian)))
recov.ci	<-	exp(matrix(c(recov.hat$par-1.96*recov.se,recov.hat$par+1.96*recov.se),ncol=3,nrow=2,byrow=T))

c1.hat	<-	exp(recov.hat$par[1])
c2.hat	<-	exp(recov.hat$par[2])
c3.hat	<-	exp(recov.hat$par[3])

# Estimación del número de recuperados
rec.t	<-	recov.t(c1.hat,c2.hat,c3.hat,1:ncol(dat))

tmp<-pop-sims1[,"I"]-sims1[,"D"]-sims1[,"R"]

N	<-	pop
plot(dat[3,]+dat[2,],pch=19)
points(rec.t,type="l",col="red")


r0	<-	1+exp(2*c2.hat*tanh(c3.hat))
S.t	<-	N*(1-(r0/N)*rec.t+0.5*(r0/N)^2*rec.t^2)

R.t	<-	r0*S.t/N
R.t1<-	R0*sims1[,"S"]/N

I.t	<-	N-S.t
U.t	<-	I.t-rec.t
I.t1<-	I.t-U.t
prop.ind.undoc	<-	mean(U.t/I.t)


####################################################################################################
#
# Modelamiento de COVID 19 para Barranquilla y sus localidades utilizando 
# utilizando el logistico generalizado estocástico para modelar y predecir 
# el número de infectados y muertos a corto plazo
# Este archivo guarda los resultados en una lista que se llama resultados
# la cual contiene la siguiente estructura
# Seis niveles primarios para las 6 ciudades en para las que se realizó la estimación
# BARRANQUILLA, CALI, CARTAGENA, BUCARAMANGA, BOGOTA, MEDELLIN.
# Para Barranquilla, la lista contiene 7 niveles más que las demas pues tiene la estimación
# de cada una de las localidades "METROPOLITANA","NORTE - CENTRO HISTORICO","RIOMAR","SUROCCIDENTE"
# ,"SURORIENTE","NA", "Consolidado".  Para cada una de las localidades y ciudades se guardo en una 
# lista con los resultados de la estimación, predicción, proyección, riesgo y los datos 
# originales el cual contiene tanto el acumulado y el número de casos diarios. Para el nivel consolidado
# la lista contiene 3 niveles adicionales para Infectados, Fallecidos y Recuperados.
# Ver la descripción de la salida para cada una de las funciones utilizadas en los archivos del código fuente.
#
####################################################################################################

####################################################################################################
#
# Ultima actualización: 9 de Diciembre de 2020
#
####################################################################################################

####################################################################################################
# Parametros requeridos para correr el codigo
# Cargar codigo base
source(paste(path,'COVID_Baq_librerias.R',sep=''), chdir = TRUE)

# Parametros de corrida
# Descomentar esta sección y definir los parametros
#path.dat		<-	'COVID_package_v1/' #definir el directorio para guardar o buscar los datos. 
#download		<-	T #Definir si se quiere descargar los datos o buscar en el directorio.
#plot.it			<-	T	# Para generar las figuras de la estimación
#path.figs		<-	'../Figures/' # Ruta para guardar las figuras
#save.it			<-	TRUE # Para guardar los resultados compilados en la lista
#save.path		<-	'../RData/' # Ruta para guardar los resultados

resultados	<-	list() # Guarda los resultados

#######################################################################
#
# Predicción de casos infectados, Muertos y Recuperados para la ciudad
# de Barranquilla
#
#######################################################################

# Datos por INS
data.ciudad	<-	getCOVID.data(wd=path.dat
							,location="BARRANQUILLA",type="Municipio",download=download,date.type="Diagnostico"
							,data.origin="INS")

#calculo Rt
Rt_nuevos.casos		<-	data.ciudad$BARRANQUILLA[data.ciudad$BARRANQUILLA[,1]>0,2]
break.pt		<-	which(names(Rt_nuevos.casos)=="2020-09-01")	
Rt.baq		<-	rt.calc(x=Rt_nuevos.casos[1:break.pt],type="logistic")
Rt.baq2		<-	rt.calc(x=Rt_nuevos.casos[(break.pt-19):length(Rt_nuevos.casos)],type="linear")
Rt.baq		<-	rbind(Rt.baq,Rt.baq2)
Rt_fechas		<-	as.Date(names(Rt_nuevos.casos))
Rt_n.obs	<-	length(Rt_nuevos.casos)

#################################################################							
# Estimación del número de infectados fase 1.
inf.dat		<-	data.ciudad$BARRANQUILLA[data.ciudad$BARRANQUILLA[,1]>0,1]
break.pt		<-	which(names(inf.dat)=="2020-09-01")
baq.inf.mles<-	gen.log.optim(x=inf.dat[1:break.pt])
K	<-	baq.inf.mles$mles["K","MLE"]
r	<-	baq.inf.mles$mles["r","MLE"]
theta	<-	baq.inf.mles$mles["theta","MLE"]
alpha	<-	baq.inf.mles$mles["alpha","MLE"]

tot.len	<-	length(inf.dat)

obs.len		<-	length(inf.dat[1:break.pt])
baq.pred	<-	gen.log.mod(K,r,theta,alpha,t=1:obs.len)
baq.pred.ci	<-	gen.log.cis(gen.log.optim.out=baq.inf.mles,len=obs.len,Bsims=10000,sim.dist="NegBin",mean.width=20)

# Estimación del número de infectados fase 2 a partir de Septiembre 1 de 2020

new.cases			<-	data.ciudad$BARRANQUILLA[data.ciudad$BARRANQUILLA[,1]>0,2]

# Estimación del número de infectados fase 2.
# Modelando los casos nuevos
# Casos nuevos

new.cases2	<-	new.cases[(break.pt+1):(tot.len-4)]

# Modelo Poisson
x.vals		<-	1:length(new.cases2)
mod2			<-	glm(new.cases2~x.vals,family="poisson")
poisson.pred		<-	cumsum(c(tail(baq.pred,n=1),predict(mod2,type="response")))
rmse.poisson		<-	sqrt(mean((inf.dat[1:(tot.len-4)] - c(baq.pred,poisson.pred[-1]))^2))
# Modelo Negativo Binomial
mod2.1	<-	glm.nb(new.cases2~x.vals)
nb.pred	<-	cumsum(c(tail(baq.pred,n=1),predict(mod2.1,type="response")))
rmse.nb	<-	sqrt(mean((inf.dat[1:(tot.len-4)] - c(baq.pred,nb.pred[-1]))^2))

# Modelo logistico
inf.dat2		<-	inf.dat[(break.pt+1):(tot.len-4)]
mod3		<-	gen.log.optim(x=inf.dat2)
K2	<-	mod3$mles["K","MLE"]
r2	<-	mod3$mles["r","MLE"]
theta2	<-	mod3$mles["theta","MLE"]
alpha2	<-	mod3$mles["alpha","MLE"]
genlog.pred	<-	c(baq.pred,gen.log.mod(K2,r2,theta2,alpha2,t=1:length(inf.dat2)))
rmse.genlog2	<-	sqrt(mean((inf.dat[1:(tot.len-4)] - genlog.pred)^2))

# Selección del mejor modelo según el crecimiento
best.mod		<-	which.min(c(rmse.poisson,rmse.nb,rmse.genlog2))


if(best.mod==1){
	# Modelo de poisson es el mejor de acuerdo con la raíz cuadrada del error medio	
	pred.xvals	<-	1:(length(new.cases2)+34)
	mod2.pred<-	predict(mod2,newdata=data.frame(x.vals=pred.xvals),type="response")

	# Intervalos de confianza de las observaciones utilizando bootstrap parametrico
	Bsims	<-	1000
	n.cases	<-	length(new.cases2)
	boot.mat		<-	matrix(0,nrow=Bsims,ncol=length(pred.xvals))

	for(i in 1:Bsims){
	
		ith.dat	<-	rpois(n.cases,lambda=new.cases2)
		ith.mod	<-	glm(ith.dat~x.vals,family="poisson")
		ith.predict	<-	predict(ith.mod,newdata=data.frame(x.vals=pred.xvals),type="response")
		boot.mat[i,]<-	ith.predict
	
	}

	mod2.1.ci	<-	apply(boot.mat,2,quantile,prob=c(0.025,0.975))

	obs.len2	<-	length(new.cases2)+4
	Predicted.new.mat	<-	cbind(mod2.1.ci[1,],mod2.pred,mod2.1.ci[2,])
	cum.new.mat	<-	cbind(cumsum(Predicted.new.mat[,1])+tail(baq.pred.ci$Predicted[,1],n=1)
					,cumsum(Predicted.new.mat[,2])+tail(baq.pred.ci$Predicted[,2],n=1)
					,cumsum(Predicted.new.mat[,3])+tail(baq.pred.ci$Predicted[,3],n=1))
	boot.list	<-	list(Predicted.new = Predicted.new.mat[obs.len2:(obs.len2+29),], Sims=boot.mat[,obs.len2:(obs.len2+29)])
	
	# Estimación del riesgo
	baq.risk		<-	gen.log.risk(gen.log.proy.out=boot.list,threshold=50:200)
	
	# Actualización de los resultados dados los dos modelos
	baq.pred.ci$Predicted	<-	rbind(baq.pred.ci$Predicted,cum.new.mat[1:obs.len2,])
	baq.pred.ci$Predicted.new	<-	rbind(baq.pred.ci$Predicted.new,Predicted.new.mat[1:obs.len2,])
	baq.proy		<-	list(Predicted=cum.new.mat[(obs.len2):(obs.len2+29),],Predicted.new=Predicted.new.mat[(obs.len2):(obs.len2+29),])
	
		
}else if(best.mod==2){
	# Modelo negativo binomial es el mejor de acuerdo con la raíz cuadrada del error medio
	pred.xvals	<-	1:(length(new.cases2)+34)
	mod2.pred<-	predict(mod2.1,newdata=data.frame(x.vals=pred.xvals),type="response")

	# Intervalos de confianza de las observaciones utilizando bootstrap parametrico
	Bsims	<-	1000
	n.cases	<-	length(new.cases2)
	boot.mat		<-	matrix(0,nrow=Bsims,ncol=length(pred.xvals))

	for(i in 1:Bsims){
	
		ith.dat	<-	rnbinom(n.cases,mu=new.cases2,size=mod2.1$theta)
		ith.mod	<-	glm.nb(ith.dat~x.vals)
		ith.predict	<-	predict(ith.mod,newdata=data.frame(x.vals=pred.xvals),type="response")
		boot.mat[i,]<-	ith.predict
	
	}

	mod2.1.ci	<-	apply(boot.mat,2,quantile,prob=c(0.025,0.975))

	obs.len2	<-	length(new.cases2)+4
	Predicted.new.mat	<-	cbind(mod2.1.ci[1,],mod2.pred,mod2.1.ci[2,])
	cum.new.mat	<-	cbind(cumsum(Predicted.new.mat[,1])+tail(baq.pred.ci$Predicted[,1],n=1)
						,cumsum(Predicted.new.mat[,2])+tail(baq.pred.ci$Predicted[,2],n=1)
						,cumsum(Predicted.new.mat[,3])+tail(baq.pred.ci$Predicted[,3],n=1))
	boot.list	<-	list(Predicted.new = Predicted.new.mat[obs.len2:(obs.len2+29),], Sims=boot.mat[,obs.len2:(obs.len2+29)])
	
	# Estimación del riesgo
	baq.risk		<-	gen.log.risk(gen.log.proy.out=boot.list,threshold=50:200)
	
	# Actualización de los resultados dados los dos modelos
	baq.pred.ci$Predicted	<-	rbind(baq.pred.ci$Predicted,cum.new.mat[1:obs.len2,])
	baq.pred.ci$Predicted.new	<-	rbind(baq.pred.ci$Predicted.new,Predicted.new.mat[1:obs.len2,])
	baq.proy		<-	list(Predicted=cum.new.mat[(obs.len2):(obs.len2+29),],Predicted.new=Predicted.new.mat[(obs.len2):(obs.len2+29),])
	
} else{
	# Modelo logistico generalizado es el mejor de acuerdo con la raíz cuadrada del error medio
	pred.xvals	<-	1:(length(new.cases2)+4)
	mod3.pred	<-	gen.log.mod(K2,r2,theta2,alpha2,t=pred.xvals)
	new.cases.init	<-	new.cases2[1]
	mod3.ci		<-	tryCatch(gen.log.cis(gen.log.optim.out=mod3,len=length(pred.xvals),Bsims=10000,sim.dist="NegBin",mean.width=20,n0=new.cases.init)
							,error=function(e){gen.log.cis(gen.log.optim.out=mod3,len=length(pred.xvals),Bsims=10000,sim.dist="Poisson",n0=new.cases.init)})
	obs.len2	<-	length(new.cases2)+4
	mod3.proy	<-	tryCatch(gen.log.predict(mod3,t0=obs.len2,len=30,Bsims=10000,sim.dist="NegBin",mean.width=20,n0=new.cases.init)
							,error=function(e)gen.log.predict(mod3,t0=obs.len2,len=30,Bsims=10000,sim.dist="Poisson",n0=new.cases.init))
							
	# Actualización de los resultados dados los dos modelos
	baq.pred.ci$Predicted	<-	rbind(baq.pred.ci$Predicted,mod3.ci$Predicted)
	baq.pred.ci$Predicted.new	<-	rbind(baq.pred.ci$Predicted.new,mod3.ci$Predicted.new)
	baq.proy		<-	mod3.proy
	
	# Estimación del riesgo
	baq.risk		<-	gen.log.risk(gen.log.proy.out=baq.proy,threshold=50:200)
	
}

first.date	<-	head(names(inf.dat),n=1)
obs.dates	<-	seq.Date(as.Date(first.date),by="days",length.out=tot.len)
last.date	<-	tail(names(inf.dat),n=1)
proy.dates<-seq.Date(as.Date(last.date),by="days",length.out=30)


Datos=data.frame(acumulado=inf.dat,dario=new.cases)

resultados$Barranquilla$Consolidado$Infectados	<-	list(Estimacion=baq.inf.mles,Prediccion=baq.pred.ci
														,Proyeccion=baq.proy
														,Riesgo=baq.risk,fechas=list(obs=obs.dates,proy=proy.dates),Datos=Datos)

resultados$Barranquilla$Consolidado$Rt <- list(Rt=Rt.baq,fechas=Rt_fechas,nobs=Rt_n.obs, nuevos_casos=Rt_nuevos.casos)

#################################################################
# Figuras de la estimación
first.date	<-	head(names(inf.dat),n=1)
obs.dates	<-	seq.Date(as.Date(first.date),by="days",length.out=tot.len)
last.date	<-	tail(names(inf.dat),n=1)
proy.dates<-seq.Date(as.Date(last.date),by="days",length.out=30)

if(plot.it){
	fig.name		<-	paste(path.figs,"Infectados.tiff",sep="")
	tiff(fig.name,width=6,height=7,units="in",res=300
		,type="cairo",compression="lzw",family="times")
	mat<-matrix(c(1,1,2,2,
				 0,0,0,0,
				 0,3,3,0,
				 0,4,4,0),nrow=4,ncol=4,byrow=TRUE)
	layout(mat,heights=c(1,0.05,1,0.15))
	par(mar=c(2,2.5,1,1),oma=c(1,0,2,0),mgp=c(1.5,0.5,0),tcl=-0.25)
	
	# Figura A. Numero de casos nuevos observados, predichos y proyectados
	plot(obs.dates,new.cases,pch=19,xlim=range(c(obs.dates,proy.dates)),bty="l",xlab="",ylab="Número de casos nuevos"
		,cex=0.5,cex.axis=0.7,cex.lab=0.7)
	points(obs.dates,baq.pred.ci$Predicted.new[,2],type="l",col="red")
	points(obs.dates,baq.pred.ci$Predicted.new[,1],type="l",col="red",lty=2)
	points(obs.dates,baq.pred.ci$Predicted.new[,3],type="l",col="red",lty=2)
	points(proy.dates,baq.proy$Predicted.new[,1],type="l",col="darkorange",lty=2)
	points(proy.dates,baq.proy$Predicted.new[,3],type="l",col="darkorange",lty=2)
	points(proy.dates,baq.proy$Predicted.new[,2],type="l",col="darkorange")
	mtext("A.",side=3,adj=0.1)

	# Figura B. Numero de casos acumulados observados, predichos y proyectados
	plot(obs.dates,inf.dat,ylim=c(0,max(baq.proy$Predicted)),pch=16,xlim=range(c(obs.dates,proy.dates)),bty="l",xlab=""
		,ylab="Número acumulado de casos"	,cex=0.5,cex.axis=0.7,cex.lab=0.7)
	points(obs.dates,baq.pred.ci$Predicted[,2],type="l",col="red")
	points(obs.dates,baq.pred.ci$Predicted[,1],type="l",col="red",lty=2)
	points(obs.dates,baq.pred.ci$Predicted[,3],type="l",col="red",lty=2)
	points(proy.dates,baq.proy$Predicted[,1],type="l",col="darkorange",lty=2)
	points(proy.dates,baq.proy$Predicted[,2],type="l",col="darkorange")
	points(proy.dates,baq.proy$Predicted[,3],type="l",col="darkorange",lty=2)
	mtext("B.",side=3,adj=0.1)
	
	# Figura C. Probabilidad de observar un número determinado de casos nuevos
	image(x=proy.dates,y=50:200,t(baq.risk),zlim=c(0,1),xlim=range(proy.dates)
		,ylim=c(50,200),col=rev(my.col.ramp(151)),xlab="",ylab="Número de casos nuevos"
		,cex.axis=0.7,cex.lab=0.7)
	box()
	mtext("C.",side=3,adj=0.1)
	par(mar=c(2.5,2.5,0,1))
	image(x=50:200,y=1,z=matrix(seq(0,1,length.out=151),ncol=1),col=rev(my.col.ramp(151))
			,xaxt="n",yaxt="n",xlab="Probabilidad de nuevos casos",cex.lab=0.7,ylab="")
	box()
	axis(1,at=seq(50,200,length.out=5),labels=seq(0,1,length.out=5),cex.axis=0.7,cex.lab=0.7)

	# Leyenda de la figura
	par(fig=c(0,1,0,1),new=TRUE,mar=c(0,3,0,1),oma=c(0,0,0,0))
	plot(0,0,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
	legend("top",legend=c("Observado","Predicho","Proyectado","Incertidumbre en la Predicción"
						,"Incertidumbre en la Proyección")
		,pch=c(16,-1,-1,-1,-1),lty=c(-1,1,1,2,2),col=c("black","red","darkorange","red","darkorange")
		,bty="n",ncol=3,cex=0.7)
dev.off()
}

###########################################################################################
# Estimación de las muertes a través de la epidemia
mort.dat			<-	data.ciudad$BARRANQUILLA[data.ciudad$BARRANQUILLA[,3]>0,3]
break.pt			<-	which(names(mort.dat)=="2020-09-01")
tot.len			<-	length(mort.dat)
baq.mort.mles	<-	gen.log.optim(mort.dat[1:break.pt],hessian=TRUE)

K		<-	baq.mort.mles$mles["K","MLE"]
r		<-	baq.mort.mles$mles["r","MLE"]
theta	<-	baq.mort.mles$mles["theta","MLE"]
alpha	<-	baq.mort.mles$mles["alpha","MLE"]

obs.len		<-	length(mort.dat[1:break.pt])
baq.pred	<-	gen.log.mod(K,r,theta,alpha,t=1:obs.len)
baq.pred.ci	<-	gen.log.cis(baq.mort.mles,len=obs.len,sim.dist="Poisson",Bsims=10000,mean.width=20)

# Estimación del número de infectados fase 2 a partir de Septiembre 1 de 2020

new.cases			<-	data.ciudad$BARRANQUILLA[data.ciudad$BARRANQUILLA[,3]>0,4]

# Estimación del número de infectados fase 2.
# Modelando los casos nuevos
# Casos nuevos

new.cases2	<-	new.cases[(break.pt+1):(tot.len)]

# Modelo Poisson
x.vals		<-	1:length(new.cases2)
mod2			<-	glm(new.cases2~x.vals,family="poisson")
poisson.pred		<-	cumsum(c(tail(baq.pred,n=1),predict(mod2,type="response")))
rmse.poisson		<-	sqrt(mean((mort.dat[1:(tot.len)] - c(baq.pred,poisson.pred[-1]))^2))
# Modelo Negativo Binomial
mod2.1	<-	glm.nb(new.cases2~x.vals)
nb.pred	<-	cumsum(c(tail(baq.pred,n=1),predict(mod2.1,type="response")))
rmse.nb	<-	sqrt(mean((mort.dat[1:(tot.len)] - c(baq.pred,nb.pred[-1]))^2))

# Modelo logistico
mort.dat2		<-	mort.dat[(break.pt+1):(tot.len)]
mod3		<-	gen.log.optim(x=mort.dat2)
K2	<-	mod3$mles["K","MLE"]
r2	<-	mod3$mles["r","MLE"]
theta2	<-	mod3$mles["theta","MLE"]
alpha2	<-	mod3$mles["alpha","MLE"]
genlog.pred	<-	c(baq.pred,gen.log.mod(K2,r2,theta2,alpha2,t=1:length(mort.dat2)))
rmse.genlog2	<-	sqrt(mean((mort.dat[1:(tot.len)] - genlog.pred)^2))

# Selección del mejor modelo según el crecimiento
best.mod		<-	which.min(c(rmse.poisson,rmse.nb,rmse.genlog2))


if(best.mod==1){
	# Modelo de poisson es el mejor de acuerdo con la raíz cuadrada del error medio	
	pred.xvals	<-	1:(length(new.cases2)+30)
	mod2.pred<-	predict(mod2,newdata=data.frame(x.vals=pred.xvals),type="response")

	# Intervalos de confianza de las observaciones utilizando bootstrap parametrico
	Bsims	<-	1000
	n.cases	<-	length(new.cases2)
	boot.mat		<-	matrix(0,nrow=Bsims,ncol=length(pred.xvals))

	for(i in 1:Bsims){
	
		ith.dat	<-	rpois(n.cases,lambda=new.cases2)
		ith.mod	<-	glm(ith.dat~x.vals,family="poisson")
		ith.predict	<-	predict(ith.mod,newdata=data.frame(x.vals=pred.xvals),type="response")
		boot.mat[i,]<-	ith.predict
	
	}

	mod2.1.ci	<-	apply(boot.mat,2,quantile,prob=c(0.025,0.975))

	obs.len2	<-	length(new.cases2)
	Predicted.new.mat	<-	cbind(mod2.1.ci[1,],mod2.pred,mod2.1.ci[2,])
	cum.new.mat	<-	cbind(cumsum(Predicted.new.mat[,1])+tail(baq.pred.ci$Predicted[,1],n=1)
					,cumsum(Predicted.new.mat[,2])+tail(baq.pred.ci$Predicted[,2],n=1)
					,cumsum(Predicted.new.mat[,3])+tail(baq.pred.ci$Predicted[,3],n=1))
	boot.list	<-	list(Predicted.new = Predicted.new.mat[obs.len2:(obs.len2+29),], Sims=boot.mat[,obs.len2:(obs.len2+29)])
	
	# Estimación del riesgo
	baq.risk		<-	gen.log.risk(gen.log.proy.out=boot.list,threshold=1:10)
	
	# Actualización de los resultados dados los dos modelos
	baq.pred.ci$Predicted	<-	rbind(baq.pred.ci$Predicted,cum.new.mat[1:obs.len2,])
	baq.pred.ci$Predicted.new	<-	rbind(baq.pred.ci$Predicted.new,Predicted.new.mat[1:obs.len2,])
	baq.proy		<-	list(Predicted=cum.new.mat[(obs.len2):(obs.len2+29),],Predicted.new=Predicted.new.mat[(obs.len2):(obs.len2+29),])
	
		
}else if(best.mod==2){
	# Modelo negativo binomial es el mejor de acuerdo con la raíz cuadrada del error medio
	pred.xvals	<-	1:(length(new.cases2)+30)
	mod2.pred<-	predict(mod2.1,newdata=data.frame(x.vals=pred.xvals),type="response")

	# Intervalos de confianza de las observaciones utilizando bootstrap parametrico
	Bsims	<-	1000
	n.cases	<-	length(new.cases2)
	boot.mat		<-	matrix(0,nrow=Bsims,ncol=length(pred.xvals))

	for(i in 1:Bsims){
	
		ith.dat	<-	rnbinom(n.cases,mu=new.cases2,size=mod2.1$theta)
		ith.mod	<-	glm.nb(ith.dat~x.vals)
		ith.predict	<-	predict(ith.mod,newdata=data.frame(x.vals=pred.xvals),type="response")
		boot.mat[i,]<-	ith.predict
	
	}

	mod2.1.ci	<-	apply(boot.mat,2,quantile,prob=c(0.025,0.975))

	obs.len2	<-	length(new.cases2)
	Predicted.new.mat	<-	cbind(mod2.1.ci[1,],mod2.pred,mod2.1.ci[2,])
	cum.new.mat	<-	cbind(cumsum(Predicted.new.mat[,1])+tail(baq.pred.ci$Predicted[,1],n=1)
						,cumsum(Predicted.new.mat[,2])+tail(baq.pred.ci$Predicted[,2],n=1)
						,cumsum(Predicted.new.mat[,3])+tail(baq.pred.ci$Predicted[,3],n=1))
	boot.list	<-	list(Predicted.new = Predicted.new.mat[obs.len2:(obs.len2+29),], Sims=boot.mat[,obs.len2:(obs.len2+29)])
	
	# Estimación del riesgo
	baq.risk		<-	gen.log.risk(gen.log.proy.out=boot.list,threshold=1:10)
	
	# Actualización de los resultados dados los dos modelos
	baq.pred.ci$Predicted	<-	rbind(baq.pred.ci$Predicted,cum.new.mat[1:obs.len2,])
	baq.pred.ci$Predicted.new	<-	rbind(baq.pred.ci$Predicted.new,Predicted.new.mat[1:obs.len2,])
	baq.proy		<-	list(Predicted=cum.new.mat[(obs.len2):(obs.len2+29),],Predicted.new=Predicted.new.mat[(obs.len2):(obs.len2+29),])
	
} else{
	# Modelo logistico generalizado es el mejor de acuerdo con la raíz cuadrada del error medio
	pred.xvals	<-	1:(length(new.cases2))
	mod3.pred	<-	gen.log.mod(K2,r2,theta2,alpha2,t=pred.xvals)
	new.cases.init	<-	new.cases2[1]
	mod3.ci		<-	tryCatch(gen.log.cis(gen.log.optim.out=mod3,len=length(pred.xvals),Bsims=10000,sim.dist="NegBin",mean.width=20,n0=new.cases.init)
							,error=function(e){gen.log.cis(gen.log.optim.out=mod3,len=length(pred.xvals),Bsims=10000,sim.dist="Poisson",n0=new.cases.init)})
	obs.len2	<-	length(new.cases2)
	mod3.proy	<-	tryCatch(gen.log.predict(mod3,t0=obs.len2,len=30,Bsims=10000,sim.dist="NegBin",mean.width=20,n0=new.cases.init)
							,error=function(e)gen.log.predict(mod3,t0=obs.len2,len=30,Bsims=10000,sim.dist="Poisson",n0=new.cases.init))
							
	# Actualización de los resultados dados los dos modelos
	baq.pred.ci$Predicted	<-	rbind(baq.pred.ci$Predicted,mod3.ci$Predicted)
	baq.pred.ci$Predicted.new	<-	rbind(baq.pred.ci$Predicted.new,mod3.ci$Predicted.new)
	baq.proy		<-	mod3.proy
	
	# Estimación del riesgo
	baq.risk		<-	gen.log.risk(gen.log.proy.out=baq.proy,threshold=1:10)
	
}

first.date	<-	head(names(mort.dat),n=1)
obs.dates	<-	seq.Date(as.Date(first.date),by="days",length.out=tot.len)
last.date	<-	tail(names(mort.dat),n=1)
proy.dates<-seq.Date(as.Date(last.date),by="days",length.out=30)

Datos=data.frame(acumulado=mort.dat,dario=new.cases)

resultados$Barranquilla$Consolidado$Fallecidos	<-	list(Estimacion=baq.mort.mles,Prediccion=baq.pred.ci
														,Proyeccion=baq.proy
														,Riesgo=baq.risk,fechas=list(obs=obs.dates,proy=proy.dates),Datos=Datos)

###########################################################################################
# Figuras de la estimación
first.date	<-	head(names(mort.dat),n=1)
obs.dates	<-	seq.Date(as.Date(first.date),by="days",length.out=tot.len)
last.date	<-	tail(names(mort.dat),n=1)
proy.dates<-seq.Date(as.Date(last.date),by="days",length.out=30)

if(plot.it){
	fig.name		<-	paste(path.figs,"Muertes.tiff",sep="")
	tiff(fig.name,width=6,height=7,units="in",res=300
		,type="cairo",compression="lzw",family="times")
	mat<-matrix(c(1,1,2,2,
				 0,0,0,0,
				 0,3,3,0,
				 0,4,4,0),nrow=4,ncol=4,byrow=TRUE)
	layout(mat,heights=c(1,0.05,1,0.15))
	par(mar=c(2,2.5,1,1),oma=c(1,0,2,0),mgp=c(1.5,0.5,0),tcl=-0.25)
	
	# Figura A. Numero de casos nuevos observados, predichos y proyectados
	plot(obs.dates,new.cases,pch=19,xlim=range(c(obs.dates,proy.dates)),bty="l",xlab="",ylab="Número de casos nuevos"
		,cex=0.5,cex.axis=0.7,cex.lab=0.7)
	points(obs.dates,baq.pred.ci$Predicted.new[,2],type="l",col="red")
	points(obs.dates,baq.pred.ci$Predicted.new[,1],type="l",col="red",lty=2)
	points(obs.dates,baq.pred.ci$Predicted.new[,3],type="l",col="red",lty=2)
	points(proy.dates,baq.proy$Predicted.new[,1],type="l",col="darkorange",lty=2)
	points(proy.dates,baq.proy$Predicted.new[,3],type="l",col="darkorange",lty=2)
	points(proy.dates,baq.proy$Predicted.new[,2],type="l",col="darkorange")
	mtext("A.",side=3,adj=0.1)

	# Figura B. Numero de casos acumulados observados, predichos y proyectados
	plot(obs.dates,mort.dat,ylim=c(0,max(baq.proy$Predicted)),pch=16,xlim=range(c(obs.dates,proy.dates)),bty="l",xlab=""
		,ylab="Número acumulado de casos"	,cex=0.5,cex.axis=0.7,cex.lab=0.7)
	points(obs.dates,baq.pred.ci$Predicted[,2],type="l",col="red")
	points(obs.dates,baq.pred.ci$Predicted[,1],type="l",col="red",lty=2)
	points(obs.dates,baq.pred.ci$Predicted[,3],type="l",col="red",lty=2)
	points(proy.dates,baq.proy$Predicted[,1],type="l",col="darkorange",lty=2)
	points(proy.dates,baq.proy$Predicted[,2],type="l",col="darkorange")
	points(proy.dates,baq.proy$Predicted[,3],type="l",col="darkorange",lty=2)
	mtext("B.",side=3,adj=0.1)

	# Figura C. Probabilidad de observar un número determinado de casos nuevos
	image(x=proy.dates,y=1:10,t(baq.risk),zlim=c(0,1),xlim=range(proy.dates)
		,ylim=c(1,10),col=rev(my.col.ramp(100)),xlab="",ylab="Número de casos nuevos"
		,cex.axis=0.7,cex.lab=0.7)
	box()
	mtext("C.",side=3,adj=0.1)
	par(mar=c(2.5,2.5,0,1))
	image(x=1:100,y=1,z=matrix(seq(0,1,length.out=100),ncol=1),col=rev(my.col.ramp(100))
			,xaxt="n",yaxt="n",xlab="Probabilidad de nuevos casos",cex.lab=0.7,ylab="")
	box()
	axis(1,at=seq(1,100,length.out=5),labels=seq(0,1,length.out=5),cex.axis=0.7,cex.lab=0.7)

	# Leyenda de la figura
	par(fig=c(0,1,0,1),new=TRUE,mar=c(0,3,0,1),oma=c(0,0,0,0))
	plot(0,0,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
	legend("top",legend=c("Observado","Predicho","Proyectado","Incertidumbre en la Predicción"
							,"Incertidumbre en la Proyección")
						,pch=c(16,-1,-1,-1,-1),lty=c(-1,1,1,2,2)
						,col=c("black","red","darkorange","red","darkorange")
						,bty="n",ncol=3,cex=0.7)
	dev.off()
}

###########################################################################################
# Estimación del Recuperados
recov.dat		<-	data.ciudad$BARRANQUILLA[data.ciudad$BARRANQUILLA[,5]>0,5]
break.pt			<-	which(names(recov.dat)=="2020-09-01")
tot.len			<-	length(recov.dat)

baq.recov.mles	<-	gen.log.optim(recov.dat[1:break.pt],hessian=TRUE)

K	<-	baq.recov.mles$mles["K","MLE"]
r	<-	baq.recov.mles$mles["r","MLE"]
theta	<-	baq.recov.mles$mles["theta","MLE"]
alpha	<-	baq.recov.mles$mles["alpha","MLE"]

obs.len		<-	length(recov.dat[1:break.pt])
baq.pred	<-	gen.log.mod(K,r,theta,alpha,t=1:obs.len)
baq.pred.ci	<-	gen.log.cis(baq.recov.mles,len=obs.len,Bsims=10000,sim.dist="NegBin")

# Estimación del número de infectados fase 2 a partir de Septiembre 1 de 2020

new.cases			<-	data.ciudad$BARRANQUILLA[data.ciudad$BARRANQUILLA[,5]>0,6]

# Estimación del número de infectados fase 2.
# Modelando los casos nuevos
# Casos nuevos

new.cases2	<-	new.cases[(break.pt+1):(tot.len)]

# Modelo Poisson
x.vals		<-	1:length(new.cases2)
mod2			<-	glm(new.cases2~x.vals,family="poisson")
poisson.pred		<-	cumsum(c(tail(baq.pred,n=1),predict(mod2,type="response")))
rmse.poisson		<-	sqrt(mean((recov.dat[1:(tot.len)] - c(baq.pred,poisson.pred[-1]))^2))
# Modelo Negativo Binomial
mod2.1	<-	glm.nb(new.cases2~x.vals)
nb.pred	<-	cumsum(c(tail(baq.pred,n=1),predict(mod2.1,type="response")))
rmse.nb	<-	sqrt(mean((recov.dat[1:(tot.len)] - c(baq.pred,nb.pred[-1]))^2))

# Modelo logistico
recov.dat2		<-	recov.dat[(break.pt+1):(tot.len)]
mod3		<-	gen.log.optim(x=recov.dat2)
K2	<-	mod3$mles["K","MLE"]
r2	<-	mod3$mles["r","MLE"]
theta2	<-	mod3$mles["theta","MLE"]
alpha2	<-	mod3$mles["alpha","MLE"]
genlog.pred	<-	c(baq.pred,gen.log.mod(K2,r2,theta2,alpha2,t=1:length(recov.dat2)))
rmse.genlog2	<-	sqrt(mean((recov.dat[1:(tot.len)] - genlog.pred)^2))

# Selección del mejor modelo según el crecimiento
best.mod		<-	which.min(c(rmse.poisson,rmse.nb,rmse.genlog2))


if(best.mod==1){
	# Modelo de poisson es el mejor de acuerdo con la raíz cuadrada del error medio	
	pred.xvals	<-	1:(length(new.cases2)+30)
	mod2.pred<-	predict(mod2,newdata=data.frame(x.vals=pred.xvals),type="response")

	# Intervalos de confianza de las observaciones utilizando bootstrap parametrico
	Bsims	<-	1000
	n.cases	<-	length(new.cases2)
	boot.mat		<-	matrix(0,nrow=Bsims,ncol=length(pred.xvals))

	for(i in 1:Bsims){
	
		ith.dat	<-	rpois(n.cases,lambda=new.cases2)
		ith.mod	<-	glm(ith.dat~x.vals,family="poisson")
		ith.predict	<-	predict(ith.mod,newdata=data.frame(x.vals=pred.xvals),type="response")
		boot.mat[i,]<-	ith.predict
	
	}

	mod2.1.ci	<-	apply(boot.mat,2,quantile,prob=c(0.025,0.975))

	obs.len2	<-	length(new.cases2)
	Predicted.new.mat	<-	cbind(mod2.1.ci[1,],mod2.pred,mod2.1.ci[2,])
	cum.new.mat	<-	cbind(cumsum(Predicted.new.mat[,1])+tail(baq.pred.ci$Predicted[,1],n=1)
					,cumsum(Predicted.new.mat[,2])+tail(baq.pred.ci$Predicted[,2],n=1)
					,cumsum(Predicted.new.mat[,3])+tail(baq.pred.ci$Predicted[,3],n=1))
	boot.list	<-	list(Predicted.new = Predicted.new.mat[obs.len2:(obs.len2+29),], Sims=boot.mat[,obs.len2:(obs.len2+29)])
	
	# Estimación del riesgo
	baq.risk		<-	gen.log.risk(gen.log.proy.out=boot.list,threshold=1:100)
	
	# Actualización de los resultados dados los dos modelos
	baq.pred.ci$Predicted	<-	rbind(baq.pred.ci$Predicted,cum.new.mat[1:obs.len2,])
	baq.pred.ci$Predicted.new	<-	rbind(baq.pred.ci$Predicted.new,Predicted.new.mat[1:obs.len2,])
	baq.proy		<-	list(Predicted=cum.new.mat[(obs.len2):(obs.len2+29),],Predicted.new=Predicted.new.mat[(obs.len2):(obs.len2+29),])
	
		
}else if(best.mod==2){
	# Modelo negativo binomial es el mejor de acuerdo con la raíz cuadrada del error medio
	pred.xvals	<-	1:(length(new.cases2)+30)
	mod2.pred<-	predict(mod2.1,newdata=data.frame(x.vals=pred.xvals),type="response")

	# Intervalos de confianza de las observaciones utilizando bootstrap parametrico
	Bsims	<-	1000
	n.cases	<-	length(new.cases2)
	boot.mat		<-	matrix(0,nrow=Bsims,ncol=length(pred.xvals))

	for(i in 1:Bsims){
	
		ith.dat	<-	rnbinom(n.cases,mu=new.cases2,size=mod2.1$theta)
		ith.mod	<-	glm.nb(ith.dat~x.vals)
		ith.predict	<-	predict(ith.mod,newdata=data.frame(x.vals=pred.xvals),type="response")
		boot.mat[i,]<-	ith.predict
	
	}

	mod2.1.ci	<-	apply(boot.mat,2,quantile,prob=c(0.025,0.975))

	obs.len2	<-	length(new.cases2)
	Predicted.new.mat	<-	cbind(mod2.1.ci[1,],mod2.pred,mod2.1.ci[2,])
	cum.new.mat	<-	cbind(cumsum(Predicted.new.mat[,1])+tail(baq.pred.ci$Predicted[,1],n=1)
						,cumsum(Predicted.new.mat[,2])+tail(baq.pred.ci$Predicted[,2],n=1)
						,cumsum(Predicted.new.mat[,3])+tail(baq.pred.ci$Predicted[,3],n=1))
	boot.list	<-	list(Predicted.new = Predicted.new.mat[obs.len2:(obs.len2+29),], Sims=boot.mat[,obs.len2:(obs.len2+29)])
	
	# Estimación del riesgo
	baq.risk		<-	gen.log.risk(gen.log.proy.out=boot.list,threshold=1:100)
	
	# Actualización de los resultados dados los dos modelos
	baq.pred.ci$Predicted	<-	rbind(baq.pred.ci$Predicted,cum.new.mat[1:obs.len2,])
	baq.pred.ci$Predicted.new	<-	rbind(baq.pred.ci$Predicted.new,Predicted.new.mat[1:obs.len2,])
	baq.proy		<-	list(Predicted=cum.new.mat[(obs.len2):(obs.len2+29),],Predicted.new=Predicted.new.mat[(obs.len2):(obs.len2+29),])
	
} else{
	# Modelo logistico generalizado es el mejor de acuerdo con la raíz cuadrada del error medio
	pred.xvals	<-	1:(length(new.cases2))
	mod3.pred	<-	gen.log.mod(K2,r2,theta2,alpha2,t=pred.xvals)
	new.cases.init	<-	new.cases2[1]
	mod3.ci		<-	tryCatch(gen.log.cis(gen.log.optim.out=mod3,len=length(pred.xvals),Bsims=10000,sim.dist="NegBin",mean.width=20,n0=new.cases.init)
							,error=function(e){gen.log.cis(gen.log.optim.out=mod3,len=length(pred.xvals),Bsims=10000,sim.dist="Poisson",n0=new.cases.init)})
	obs.len2	<-	length(new.cases2)
	mod3.proy	<-	tryCatch(gen.log.predict(mod3,t0=obs.len2,len=30,Bsims=10000,sim.dist="NegBin",mean.width=20,n0=new.cases.init)
							,error=function(e)gen.log.predict(mod3,t0=obs.len2,len=30,Bsims=10000,sim.dist="Poisson",n0=new.cases.init))
							
	# Actualización de los resultados dados los dos modelos
	baq.pred.ci$Predicted	<-	rbind(baq.pred.ci$Predicted,mod3.ci$Predicted)
	baq.pred.ci$Predicted.new	<-	rbind(baq.pred.ci$Predicted.new,mod3.ci$Predicted.new)
	baq.proy		<-	mod3.proy
	
	# Estimación del riesgo
	baq.risk		<-	gen.log.risk(gen.log.proy.out=baq.proy,threshold=1:100)
	
}


first.date	<-	head(names(recov.dat),n=1)
obs.dates	<-	seq.Date(as.Date(first.date),by="days",length.out=tot.len)
last.date	<-	tail(names(recov.dat),n=1)
proy.dates<-seq.Date(as.Date(last.date),by="days",length.out=30)

Datos=data.frame(acumulado=recov.dat,dario=new.cases)

resultados$Barranquilla$Consolidado$Recuperados	<-	list(Estimacion=baq.recov.mles,Prediccion=baq.pred.ci
														,Proyeccion=baq.proy
														,Riesgo=baq.risk,fechas=list(obs=obs.dates,proy=proy.dates),Datos=Datos)
														
###########################################################################################
# Figuras de la estimación
first.date	<-	head(names(recov.dat),n=1)
obs.dates	<-	seq.Date(as.Date(first.date),by="days",length.out=tot.len)
last.date	<-	tail(names(recov.dat),n=1)
proy.dates<-seq.Date(as.Date(last.date),by="days",length.out=30)

if(plot.it){
	fig.name		<-	paste(path.figs,"Recuperados.tiff",sep="")
	tiff(fig.name,width=6,height=7,units="in",res=300
		,type="cairo",compression="lzw",family="times")
	mat<-matrix(c(1,1,2,2,
				 0,0,0,0,
				 0,3,3,0,
				 0,4,4,0),nrow=4,ncol=4,byrow=TRUE)
	layout(mat,heights=c(1,0.05,1,0.15))
	par(mar=c(2,2.5,1,1),oma=c(1,0,2,0),mgp=c(1.5,0.5,0),tcl=-0.25)
	
	# Figura A. Numero de casos nuevos observados, predichos y proyectados
	plot(obs.dates,new.cases,pch=19,xlim=range(c(obs.dates,proy.dates)),bty="l",xlab="",ylab="Número de casos nuevos"
		,cex=0.5,cex.axis=0.7,cex.lab=0.7)
	points(obs.dates,baq.pred.ci$Predicted.new[,2],type="l",col="red")
	points(obs.dates,baq.pred.ci$Predicted.new[,1],type="l",col="red",lty=2)
	points(obs.dates,baq.pred.ci$Predicted.new[,3],type="l",col="red",lty=2)
	points(proy.dates,baq.proy$Predicted.new[,1],type="l",col="darkorange",lty=2)
	points(proy.dates,baq.proy$Predicted.new[,3],type="l",col="darkorange",lty=2)
	points(proy.dates,baq.proy$Predicted.new[,2],type="l",col="darkorange")
	mtext("A.",side=3,adj=0.1)

	# Figura B. Numero de casos acumulados observados, predichos y proyectados
	plot(obs.dates,recov.dat,ylim=c(0,max(baq.proy$Predicted)),pch=16,xlim=range(c(obs.dates,proy.dates))
		,bty="l",xlab=""
		,ylab="Número acumulado de casos"	,cex=0.5,cex.axis=0.7,cex.lab=0.7)
	points(obs.dates,baq.pred.ci$Predicted[,2],type="l",col="red")
	points(obs.dates,baq.pred.ci$Predicted[,1],type="l",col="red",lty=2)
	points(obs.dates,baq.pred.ci$Predicted[,3],type="l",col="red",lty=2)
	points(proy.dates,baq.proy$Predicted[,1],type="l",col="darkorange",lty=2)
	points(proy.dates,baq.proy$Predicted[,2],type="l",col="darkorange")
	points(proy.dates,baq.proy$Predicted[,3],type="l",col="darkorange",lty=2)
	mtext("B.",side=3,adj=0.1)

	# Figura C. Probabilidad de observar un número determinado de casos nuevos
	image(x=proy.dates,y=1:100,t(baq.risk),zlim=c(0,1),xlim=range(proy.dates),ylim=c(1,100)
			,col=rev(my.col.ramp(100)),xlab="",ylab="Número de casos nuevos"
			,cex.axis=0.7,cex.lab=0.7)
	box()
	mtext("C.",side=3,adj=0.1)
	par(mar=c(2.5,2.5,0,1))
	image(x=1:100,y=1,z=matrix(seq(0,1,length.out=100),ncol=1),col=rev(my.col.ramp(100))
			,xaxt="n",yaxt="n",xlab="Probabilidad de 	nuevos casos",cex.lab=0.7,ylab="")
	box()
	axis(1,at=seq(1,100,length.out=5),labels=seq(0,1,length.out=5),cex.axis=0.7,cex.lab=0.7)
	
	# Leyenda de la figura
	par(fig=c(0,1,0,1),new=TRUE,mar=c(0,3,0,1),oma=c(0,0,0,0))
	plot(0,0,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
	legend("top",legend=c("Observado","Predicho","Proyectado","Incertidumbre en la Predicción"
			,"Incertidumbre en la Proyección")
			,pch=c(16,-1,-1,-1,-1),lty=c(-1,1,1,2,2),col=c("black","red","darkorange","red","darkorange")
			,bty="n",ncol=3,cex=0.7)
dev.off()
}

#######################################################################
#
# Predicción de casos infectados por localidades
#
#######################################################################
alcaldia.locs	<-	c("METROPOLITANA","NORTE - CENTRO HISTORICO","RIOMAR","SUROCCIDENTE","SURORIENTE","NA")			
data.locs	<-	getCOVID.data(wd=paste(path.dat,"Data",sep=""),location=alcaldia.locs,type=NA
								 ,download=download,date.type="Diagnostico"
								 ,data.origin="Alcaldia")

for (i in 1:length(alcaldia.locs)){
  #calculo Rt
	inf.dat		<-	data.locs[[i]][data.locs[[i]][,1]>0,1]

	break.pt		<-	which(names(inf.dat)=="2020-09-01")	
	Rt_nuevos.casos		<-	data.locs[[i]][data.locs[[i]][,1]>0,2]
	Rt.baq		<-	rt.calc(x=Rt_nuevos.casos[1:break.pt],type="logistic")
	Rt.baq2		<-	rt.calc(x=Rt_nuevos.casos[(break.pt-19):length(inf.dat)],type="linear")
	Rt.baq		<-	rbind(Rt.baq,Rt.baq2)
	Rt_fechas		<-	as.Date(names(Rt_nuevos.casos))
	Rt_n.obs	<-	length(Rt_nuevos.casos)
	

	ith.inf.mles<-	gen.log.optim(x=inf.dat[1:break.pt])
	K	<-	ith.inf.mles$mles["K","MLE"]
	r	<-	ith.inf.mles$mles["r","MLE"]
	theta	<-	ith.inf.mles$mles["theta","MLE"]
	alpha	<-	ith.inf.mles$mles["alpha","MLE"]

	tot.len	<-	length(inf.dat)

	obs.len		<-	length(inf.dat[1:break.pt])
	ithloc.pred	<-	gen.log.mod(K,r,theta,alpha,t=1:obs.len)
	ithloc.pred.ci	<-	gen.log.cis(gen.log.optim.out=ith.inf.mles,len=obs.len,Bsims=10000,sim.dist="NegBin",mean.width=20)

	# Estimación del número de infectados fase 2.
	# Modelando los casos nuevos
	# Casos nuevos

	new.cases		<-	data.locs[[i]][data.locs[[i]][,1]>0,2]

	new.cases2	<-	new.cases[(break.pt+1):(tot.len-4)]

	# Modelo Poisson
	x.vals		<-	1:length(new.cases2)
	mod2		<-	glm(new.cases2~x.vals,family="poisson")
	poisson.pred		<-	cumsum(c(tail(ithloc.pred,n=1),predict(mod2,type="response")))
	rmse.poisson	<-	sqrt(mean((inf.dat[1:(tot.len-4)] - c(ithloc.pred,poisson.pred[-1]))^2))

	# Modelo Negativo Binomial
	mod2.1	<-	glm.nb(new.cases2~x.vals)
	nb.pred	<-	cumsum(c(tail(ithloc.pred,n=1),predict(mod2.1,type="response")))
	rmse.nb	<-	sqrt(mean((inf.dat[1:(tot.len-4)] - c(ithloc.pred,nb.pred[-1]))^2))	
	
	# Modelo logistico
	inf.dat2	<-	inf.dat[(break.pt+1):(tot.len-4)]
	mod3		<-	gen.log.optim(x=inf.dat2)
	K2		<-	mod3$mles["K","MLE"]
	r2		<-	mod3$mles["r","MLE"]
	theta2	<-	mod3$mles["theta","MLE"]
	alpha2	<-	mod3$mles["alpha","MLE"]
	genlog.pred	<-	c(ithloc.pred,gen.log.mod(K2,r2,theta2,alpha2,t=1:length(inf.dat2)))
	rmse.genlog2<-	sqrt(mean((inf.dat[1:(tot.len-4)] - genlog.pred)^2))

	# Selección del mejor modelo según el crecimiento
	best.mod		<-	which.min(c(rmse.poisson,rmse.nb,rmse.genlog2))

	if(best.mod==1){
		
		# Predicción del modelo negativo binomial 30 días en el futuro
		pred.xvals	<-	1:(length(new.cases2)+34)
		mod2.pred<-	predict(mod2,newdata=data.frame(x.vals=pred.xvals),type="response")

		# Intervalos de confianza de las observaciones utilizando bootstrap parametrico
		Bsims	<-	1000
		n.cases	<-	length(new.cases2)
		boot.mat		<-	matrix(0,nrow=Bsims,ncol=length(pred.xvals))

		for(j in 1:Bsims){
	
			ith.dat	<-	rpois(n.cases,lambda=new.cases2)
			ith.mod	<-	glm(ith.dat~x.vals,family="poisson")
			ith.predict	<-	predict(ith.mod,newdata=data.frame(x.vals=pred.xvals),type="response")
			boot.mat[j,]<-	ith.predict

		}

		mod2.1.ci	<-	apply(boot.mat,2,quantile,prob=c(0.025,0.975))

		obs.len2		<-	length(new.cases2)+4
		Predicted.new.mat	<-	cbind(mod2.1.ci[1,],mod2.pred,mod2.1.ci[2,])
		cum.new.mat	<-	cbind(cumsum(Predicted.new.mat[,1])+tail(ithloc.pred.ci$Predicted[,1],n=1)
							,cumsum(Predicted.new.mat[,2])+tail(ithloc.pred.ci$Predicted[,2],n=1)
							,cumsum(Predicted.new.mat[,3])+tail(ithloc.pred.ci$Predicted[,3],n=1))
		boot.list	<-	list(Predicted.new = Predicted.new.mat[obs.len2:(obs.len2+29),], Sims=boot.mat[,obs.len2:(obs.len2+29)])	

		# Estimación del riesgo
		ithloc.risk	<-	gen.log.risk(gen.log.proy.out=boot.list,threshold=50:300)
		
		# Actualización de los resultados dados los dos modelos
		ithloc.pred.ci$Predicted	<-	rbind(ithloc.pred.ci$Predicted,cum.new.mat[1:obs.len2,])
		ithloc.pred.ci$Predicted.new	<-	rbind(ithloc.pred.ci$Predicted.new,Predicted.new.mat[1:obs.len2,])
		ithloc.proy		<-	list(Predicted=cum.new.mat[(obs.len2):(obs.len2+29),],Predicted.new=Predicted.new.mat[(obs.len2):(obs.len2+29),])
	
	}else if(best.mod==2){
			# Predicción del modelo negativo binomial 30 días en el futuro
		pred.xvals	<-	1:(length(new.cases2)+34)
		mod2.pred<-	predict(mod2.1,newdata=data.frame(x.vals=pred.xvals),type="response")

		# Intervalos de confianza de las observaciones utilizando bootstrap parametrico	
		Bsims	<-	1000
		n.cases	<-	length(new.cases2)
		boot.mat		<-	matrix(0,nrow=Bsims,ncol=length(pred.xvals))

		for(j in 1:Bsims){
	
			ith.dat	<-	rnbinom(n.cases,mu=new.cases2,size=mod2.1$theta)
			ith.mod	<-	glm.nb(ith.dat~x.vals)
			ith.predict	<-	predict(ith.mod,newdata=data.frame(x.vals=pred.xvals),type="response")
			boot.mat[j,]<-	ith.predict

		}

		mod2.1.ci	<-	apply(boot.mat,2,quantile,prob=c(0.025,0.975))

		obs.len2		<-	length(new.cases2)+4
		Predicted.new.mat	<-	cbind(mod2.1.ci[1,],mod2.pred,mod2.1.ci[2,])
		cum.new.mat	<-	cbind(cumsum(Predicted.new.mat[,1])+tail(ithloc.pred.ci$Predicted[,1],n=1)
							,cumsum(Predicted.new.mat[,2])+tail(ithloc.pred.ci$Predicted[,2],n=1)
							,cumsum(Predicted.new.mat[,3])+tail(ithloc.pred.ci$Predicted[,3],n=1))
		boot.list	<-	list(Predicted.new = Predicted.new.mat[obs.len2:(obs.len2+29),], Sims=boot.mat[,obs.len2:(obs.len2+29)])	

		# Estimación del riesgo
		ithloc.risk	<-	gen.log.risk(gen.log.proy.out=boot.list,threshold=50:300)
		
		# Actualización de los resultados dados los dos modelos
		ithloc.pred.ci$Predicted	<-	rbind(ithloc.pred.ci$Predicted,cum.new.mat[1:obs.len2,])
		ithloc.pred.ci$Predicted.new	<-	rbind(ithloc.pred.ci$Predicted.new,Predicted.new.mat[1:obs.len2,])
		ithloc.proy		<-	list(Predicted=cum.new.mat[(obs.len2):(obs.len2+29),],Predicted.new=Predicted.new.mat[(obs.len2):(obs.len2+29),])
		

	}else{
		
		# Modelo logistico generalizado es el mejor de acuerdo con la raíz cuadrada del error medio
		pred.xvals	<-	1:(length(new.cases2)+4)
		mod3.pred	<-	gen.log.mod(K2,r2,theta2,alpha2,t=pred.xvals)
		new.cases.init	<-	new.cases2[1]
		mod3.ci		<-	tryCatch(gen.log.cis(gen.log.optim.out=mod3,len=length(pred.xvals),Bsims=10000,sim.dist="NegBin",mean.width=20,n0=new.cases.init)
								,error=function(e){gen.log.cis(gen.log.optim.out=mod3,len=length(pred.xvals),Bsims=10000,sim.dist="Poisson",n0=new.cases.init)})
		obs.len2	<-	length(new.cases2)+4
		mod3.proy	<-	tryCatch(gen.log.predict(mod3,t0=obs.len2,len=30,Bsims=10000,sim.dist="NegBin",mean.width=20,n0=new.cases.init)
								,error=function(e)gen.log.predict(mod3,t0=obs.len2,len=30,Bsims=10000,sim.dist="Poisson",n0=new.cases.init))
							
		# Actualización de los resultados dados los dos modelos
		ithloc.pred.ci$Predicted	<-	rbind(ithloc.pred.ci$Predicted,mod3.ci$Predicted)
		ithloc.pred.ci$Predicted.new	<-	rbind(ithloc.pred.ci$Predicted.new,mod3.ci$Predicted.new)
		ithloc.proy		<-	mod3.proy
	
		# Estimación del riesgo
		ithloc.risk		<-	gen.log.risk(gen.log.proy.out=ithloc.proy,threshold=50:300)
	
	}
	
	first.date	<-	head(names(inf.dat),n=1)
	obs.dates	<-	seq.Date(as.Date(first.date),by="days",length.out=tot.len)
	last.date	<-	tail(names(inf.dat),n=1)
	proy.dates<-seq.Date(as.Date(last.date),by="days",length.out=30)
	
	Datos=data.frame(acumulado=inf.dat,dario=new.cases)
	
	resultados$Barranquilla$Localidades[[i]]	<-	list(Estimacion=ith.inf.mles,Prediccion=ithloc.pred.ci
																,Proyeccion=ithloc.proy
																,Riesgo=ithloc.risk,fechas=list(obs=obs.dates,proy=proy.dates),Datos=Datos)
	#resultados$Barranquilla$Localidades[[i]]$Rt <- list(Rt=Rt.baq,fechas=Rt_fechas,nobs=Rt_n.obs, nuevos_casos=Rt_nuevos.casos)
	
}
alcaldia.locs[6]	<-	"na"
names(resultados$Barranquilla$Localidades)	<-	alcaldia.locs

# Figuras de la estimación

if(plot.it){
	fig.name		<-	paste(path.figs,"Infectados_Localidades.tiff",sep="")
	tiff(fig.name,width=7,height=10.5,units="in",res=300
		,type="cairo",compression="lzw",family="times")
	par(mfrow=c(6,2),mar=c(2,2.5,1.5,1.5),oma=c(0,0,2,0),mgp=c(1.5,0.5,0),tcl=-0.25)							

	for(i in 1:length(alcaldia.locs)){
	
	inf.dat		<-	data.locs[[i]][data.locs[[i]][,1]>0,1]
	new.cases	<-	data.locs[[i]][data.locs[[i]][,1]>0,2]
	met.pred.ci	<-	resultados$Barranquilla$Localidades[[i]]$Prediccion
	met.proy		<-	resultados$Barranquilla$Localidades[[i]]$Proyeccion
	obs.len		<-	length(inf.dat)
	
	first.date	<-	head(names(inf.dat),n=1)
	obs.dates	<-	seq.Date(as.Date(first.date),by="days",length.out=obs.len)
	last.date	<-	tail(names(inf.dat),n=1)
	proy.dates<-seq.Date(as.Date(last.date),by="days",length.out=30)
	
	# Figura A. Numero de casos nuevos observados, predichos y proyectados
	plot(obs.dates,new.cases,pch=19,xlim=range(c(obs.dates,proy.dates)),bty="l",xlab="",ylab="Número de casos nuevos"
		,cex=0.5,cex.axis=0.7,cex.lab=0.7)
	points(obs.dates,met.pred.ci$Predicted.new[,2],type="l",col="red")
	points(obs.dates,met.pred.ci$Predicted.new[,1],type="l",col="red",lty=2)
	points(obs.dates,met.pred.ci$Predicted.new[,3],type="l",col="red",lty=2)
	points(proy.dates,met.proy$Predicted.new[,1],type="l",col="darkorange",lty=2)
	points(proy.dates,met.proy$Predicted.new[,3],type="l",col="darkorange",lty=2)
	points(proy.dates,met.proy$Predicted.new[,2],type="l",col="darkorange")
	ith.loc<-alcaldia.locs[i]
	if(i==6){ith.loc<-"NO IDENTIFICADOS"}
	mtext(ith.loc,adj=0.1,cex=0.8)

	# Figura B. Numero de casos acumulados observados, predichos y proyectados
	plot(obs.dates,inf.dat,ylim=c(0,max(met.proy$Predicted)),pch=16,xlim=range(c(obs.dates,proy.dates)),bty="l",xlab=""
		,ylab="Número acumulado de casos"	,cex=0.5,cex.axis=0.7,cex.lab=0.7)
	points(obs.dates,met.pred.ci$Predicted[,2],type="l",col="red")
	points(obs.dates,met.pred.ci$Predicted[,1],type="l",col="red",lty=2)
	points(obs.dates,met.pred.ci$Predicted[,3],type="l",col="red",lty=2)
	points(proy.dates,met.proy$Predicted[,1],type="l",col="darkorange",lty=2)
	points(proy.dates,met.proy$Predicted[,2],type="l",col="darkorange")
	points(proy.dates,met.proy$Predicted[,3],type="l",col="darkorange",lty=2)

}
	# Leyenda de la figura
	par(fig=c(0,1,0,1),new=TRUE,mar=c(0,3,0,1),oma=c(0,0,0,0))
	plot(0,0,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
	legend("top",legend=c("Observado","Predicho","Proyectado","Incertidumbre en la Predicción"
							,"Incertidumbre en la Proyección")
			,pch=c(16,-1,-1,-1,-1),lty=c(-1,1,1,2,2),col=c("black","red","darkorange","red","darkorange")
			,bty="n",ncol=3,cex=0.7)
	dev.off()
}


#######################################################################
#
# Predicción de casos infectados por Ciudades: Cali, Cartagena,
# , Bucaramanga, Bogotá D.C., Medellín
#
#######################################################################
# Obtener datos
ciudades <- c("CALI","CARTAGENA","BUCARAMANGA","BOGOTA","MEDELLIN")
data.ciudad	<-	getCOVID.data(wd=path.dat
							,location=ciudades,type="Municipio",download=FALSE,date.type="Diagnostico"
							,data.origin="INS")

# Estimación del número de infectados
for(i in 2:(length(ciudades)+1)){
	inf.dat		<-	data.ciudad[[i-1]][data.ciudad[[i-1]][,1]>0,1]
	break.pt		<-	which(names(inf.dat)=="2020-09-01")
	inf.mles<-	gen.log.optim(x=inf.dat[1:break.pt])
	K	<-	inf.mles$mles["K","MLE"]
	r	<-	inf.mles$mles["r","MLE"]
	theta	<-	inf.mles$mles["theta","MLE"]
	alpha	<-	inf.mles$mles["alpha","MLE"]

	tot.len	<-	length(inf.dat)
	
	mean.width=40

	obs.len		<-	length(inf.dat[1:break.pt])
	if(inf.mles$mles["theta","97.5%"]+mean.width>obs.len){mean.width<-obs.len-inf.mles$mles["theta","97.5%"]}
	ithciud.pred	<-	gen.log.mod(K,r,theta,alpha,t=1:obs.len)
	ithciud.pred.ci	<-	gen.log.cis(gen.log.optim.out=inf.mles,len=obs.len,Bsims=10000,sim.dist="NegBin",mean.width=mean.width)

	# Estimación de la fase 2

	new.cases	<-	data.ciudad[[i-1]][data.ciudad[[i-1]][,1]>0,2]

	# Estimación del número de infectados fase 2.
	# Modelando los casos nuevos
	# Casos nuevos

	new.cases2	<-	new.cases[(break.pt+1):(tot.len-4)]

	# Modelo Poisson
	x.vals		<-	1:length(new.cases2)
	mod2		<-	glm(new.cases2~x.vals,family="poisson")
	poisson.pred		<-	cumsum(c(tail(ithciud.pred,n=1),predict(mod2,type="response")))
	rmse.poisson	<-	sqrt(mean((inf.dat[1:(tot.len-4)] - c(ithciud.pred,poisson.pred[-1]))^2))

	# Modelo Negativo Binomial
	mod2.1	<-	glm.nb(new.cases2~x.vals)
	nb.pred	<-	cumsum(c(tail(ithciud.pred,n=1),predict(mod2.1,type="response")))
	rmse.nb	<-	sqrt(mean((inf.dat[1:(tot.len-4)] - c(ithciud.pred,nb.pred[-1]))^2))	
	
	# Modelo logistico
	inf.dat2	<-	inf.dat[(break.pt+1):(tot.len-4)]
	mod3	<-	NA
	while(length(mod3)==1){
		mod3		<-	tryCatch(gen.log.optim(x=inf.dat2),error=function(e)NA)
	}
	K2		<-	mod3$mles["K","MLE"]
	r2		<-	mod3$mles["r","MLE"]
	theta2	<-	mod3$mles["theta","MLE"]
	alpha2	<-	mod3$mles["alpha","MLE"]
	genlog.pred	<-	c(ithciud.pred,gen.log.mod(K2,r2,theta2,alpha2,t=1:length(inf.dat2)))
	rmse.genlog2<-	sqrt(mean((inf.dat[1:(tot.len-4)] - genlog.pred)^2))

	# Selección del mejor modelo según el crecimiento
	best.mod		<-	which.min(c(rmse.poisson,rmse.nb,rmse.genlog2))

	if(best.mod==1){
	# Predicción del modelo negativo binomial 30 días en el futuro
	pred.xvals	<-	1:(length(new.cases2)+34)
	mod2.pred<-	predict(mod2,newdata=data.frame(x.vals=pred.xvals),type="response")

	# Intervalos de confianza de las observaciones utilizando bootstrap parametrico
	Bsims	<-	1000
	n.cases	<-	length(new.cases2)
	boot.mat		<-	matrix(0,nrow=Bsims,ncol=length(pred.xvals))

	for(j in 1:Bsims){
	
		ith.dat	<-	rpois(n.cases,lambda=new.cases2)
		ith.mod	<-	glm(ith.dat~x.vals,family="poisson")
		ith.predict	<-	predict(ith.mod,newdata=data.frame(x.vals=pred.xvals),type="response")
		boot.mat[j,]<-	ith.predict

	}

	mod2.1.ci	<-	apply(boot.mat,2,quantile,prob=c(0.025,0.975))

	obs.len2		<-	length(new.cases2)+4
	Predicted.new.mat	<-	cbind(mod2.1.ci[1,],mod2.pred,mod2.1.ci[2,])
	cum.new.mat	<-	cbind(cumsum(Predicted.new.mat[,1])+tail(ithciud.pred.ci$Predicted[,1],n=1)
						,cumsum(Predicted.new.mat[,2])+tail(ithciud.pred.ci$Predicted[,2],n=1)
						,cumsum(Predicted.new.mat[,3])+tail(ithciud.pred.ci$Predicted[,3],n=1))
	boot.list	<-	list(Predicted.new = Predicted.new.mat[obs.len2:(obs.len2+29),], Sims=boot.mat[,obs.len2:(obs.len2+29)])	

	first.date	<-	head(names(inf.dat),n=1)
	obs.dates	<-	seq.Date(as.Date(first.date),by="days",length.out=tot.len)
	last.date	<-	tail(names(inf.dat),n=1)
	proy.dates<-seq.Date(as.Date(last.date),by="days",length.out=30)

	# Actualización de los resultados dados los dos modelos
	ithciud.pred.ci$Predicted	<-	rbind(ithciud.pred.ci$Predicted,cum.new.mat[1:obs.len2,])
	ithciud.pred.ci$Predicted.new	<-	rbind(ithciud.pred.ci$Predicted.new,Predicted.new.mat[1:obs.len2,])
	ithciud.proy		<-	list(Predicted=cum.new.mat[(obs.len2):(obs.len2+29),],Predicted.new=Predicted.new.mat[(obs.len2):(obs.len2+29),])
	
	}else if(best.mod==2){
		pred.xvals	<-	1:(length(new.cases2)+34)
		mod2.pred<-	predict(mod2.1,newdata=data.frame(x.vals=pred.xvals),type="response")

		# Intervalos de confianza de las observaciones utilizando bootstrap parametrico
		Bsims	<-	1000
		n.cases	<-	length(new.cases2)
		boot.mat		<-	matrix(0,nrow=Bsims,ncol=length(pred.xvals))

		for(j in 1:Bsims){
	
			ith.dat	<-	rnbinom(n.cases,mu=new.cases2,size=mod2.1$theta)
			ith.mod	<-	glm.nb(ith.dat~x.vals)
			ith.predict	<-	predict(ith.mod,newdata=data.frame(x.vals=pred.xvals),type="response")
			boot.mat[j,]<-	ith.predict

		}

		mod2.1.ci	<-	apply(boot.mat,2,quantile,prob=c(0.025,0.975))

		obs.len2		<-	length(new.cases2)+4
		Predicted.new.mat	<-	cbind(mod2.1.ci[1,],mod2.pred,mod2.1.ci[2,])
		cum.new.mat	<-	cbind(cumsum(Predicted.new.mat[,1])+tail(ithciud.pred.ci$Predicted[,1],n=1)
							,cumsum(Predicted.new.mat[,2])+tail(ithciud.pred.ci$Predicted[,2],n=1)
							,cumsum(Predicted.new.mat[,3])+tail(ithciud.pred.ci$Predicted[,3],n=1))
		boot.list	<-	list(Predicted.new = Predicted.new.mat[obs.len2:(obs.len2+29),], Sims=boot.mat[,obs.len2:(obs.len2+29)])	

		first.date	<-	head(names(inf.dat),n=1)
		obs.dates	<-	seq.Date(as.Date(first.date),by="days",length.out=tot.len)
		last.date	<-	tail(names(inf.dat),n=1)
		proy.dates<-seq.Date(as.Date(last.date),by="days",length.out=30)

		# Actualización de los resultados dados los dos modelos
		ithciud.pred.ci$Predicted	<-	rbind(ithciud.pred.ci$Predicted,cum.new.mat[1:obs.len2,])
		ithciud.pred.ci$Predicted.new	<-	rbind(ithciud.pred.ci$Predicted.new,Predicted.new.mat[1:obs.len2,])
		ithciud.proy		<-	list(Predicted=cum.new.mat[(obs.len2):(obs.len2+29),],Predicted.new=Predicted.new.mat[(obs.len2):(obs.len2+29),])
			
	}else{
	
		# Modelo logistico generalizado es el mejor de acuerdo con la raíz cuadrada del error medio
		pred.xvals	<-	1:(length(new.cases2)+4)
		mod3.pred	<-	gen.log.mod(K2,r2,theta2,alpha2,t=pred.xvals)
		new.cases.init	<-	new.cases2[1]
		mod3.ci		<-	tryCatch(gen.log.cis(gen.log.optim.out=mod3,len=length(pred.xvals),Bsims=10000,sim.dist="NegBin"
											,mean.width=40,n0=new.cases.init)
								,error=function(e){gen.log.cis(gen.log.optim.out=mod3,len=length(pred.xvals),Bsims=10000,sim.dist="Poisson"
															,n0=new.cases.init)})
		obs.len2	<-	length(new.cases2)+4
		mod3.proy	<-	tryCatch(gen.log.predict(gen.log.optim.out=mod3,t0=obs.len2,len=30,Bsims=10000,sim.dist="NegBin",mean.width=40
												,n0=new.cases.init)
								,error=function(e)gen.log.predict(mod3,t0=obs.len2,len=30,Bsims=10000,sim.dist="Poisson"
																,n0=new.cases.init))
							
		# Actualización de los resultados dados los dos modelos
		ithciud.pred.ci$Predicted	<-	rbind(ithciud.pred.ci$Predicted,mod3.ci$Predicted)
		ithciud.pred.ci$Predicted.new	<-	rbind(ithciud.pred.ci$Predicted.new,mod3.ci$Predicted.new)
		ithciud.proy		<-	mod3.proy
	
		first.date	<-	head(names(inf.dat),n=1)
		obs.dates	<-	seq.Date(as.Date(first.date),by="days",length.out=tot.len)
		last.date	<-	tail(names(inf.dat),n=1)
		proy.dates<-seq.Date(as.Date(last.date),by="days",length.out=30)
					
	}
	
	Datos=data.frame(acumulado=inf.dat,dario=new.cases)
	
	resultados[[i]]<-	list(Estimacion=inf.mles,Prediccion=ithciud.pred.ci
							,Proyeccion=ithciud.proy,fechas=list(obs=obs.dates,proy=proy.dates),Datos=Datos)
}

names(resultados)	<-	c("BARRANQUILLA",ciudades)

if(plot.it){
	ciudades.plot	<-	c("Cali","Cartagena","Bucaramanga","Bogotá D.C.", "Medellín")
	fig.name			<-	paste(path.figs,"Infectados_Ciudades.tiff",sep="")

	tiff(fig.name,width=7,height=10.5,units="in",res=300
	,type="cairo",compression="lzw",family="times")
	par(mfrow=c(5,2),mar=c(2,2.5,1,1),oma=c(1,1,2,1),mgp=c(1.5,0.5,0),tcl=-0.25)							

	for (i in 2:(length(ciudades)+1)){
		inf.dat		<-	data.ciudad[[i-1]][data.ciudad[[i-1]][,1]>0,1]
		new.cases	<-	data.ciudad[[i-1]][data.ciudad[[i-1]][,1]>0,2]
		obs.len		<-	length(inf.dat)

		met.pred.ci	<-	resultados[[i]]$Prediccion
		met.proy		<-	resultados[[i]]$Proyeccion
	
		# Figuras de la estimación
		first.date	<-	head(names(inf.dat),n=1)
		obs.dates	<-	seq.Date(as.Date(first.date),by="days",length.out=obs.len)
		last.date	<-	tail(names(inf.dat),n=1)
		proy.dates<-seq.Date(as.Date(last.date),by="days",length.out=30)

		# Figura A. Numero de casos nuevos observados, predichos y proyectados
		plot(obs.dates,new.cases,pch=19,xlim=range(c(obs.dates,proy.dates)),bty="l",xlab="",ylab="Número de casos nuevos"
			,cex=0.5,cex.axis=0.7,cex.lab=0.7)
		points(obs.dates,met.pred.ci$Predicted.new[,2],type="l",col="red")
		points(obs.dates,met.pred.ci$Predicted.new[,1],type="l",col="red",lty=2)
		points(obs.dates,met.pred.ci$Predicted.new[,3],type="l",col="red",lty=2)
		points(proy.dates,met.proy$Predicted.new[,1],type="l",col="darkorange",lty=2)
		points(proy.dates,met.proy$Predicted.new[,3],type="l",col="darkorange",lty=2)
		points(proy.dates,met.proy$Predicted.new[,2],type="l",col="darkorange")
		ith.loc<-ciudades.plot[i-1]
		mtext(ith.loc,adj=0.1,cex=0.8)

		# Figura B. Numero de casos acumulados observados, predichos y proyectados
		plot(obs.dates,inf.dat,ylim=c(0,max(met.proy$Predicted)),pch=16,xlim=range(c(obs.dates,proy.dates)),bty="l",xlab=""
			,ylab="Número acumulado de casos"	,cex=0.5,cex.axis=0.7,cex.lab=0.7)
		points(obs.dates,met.pred.ci$Predicted[,2],type="l",col="red")
		points(obs.dates,met.pred.ci$Predicted[,1],type="l",col="red",lty=2)
		points(obs.dates,met.pred.ci$Predicted[,3],type="l",col="red",lty=2)
		points(proy.dates,met.proy$Predicted[,1],type="l",col="darkorange",lty=2)
		points(proy.dates,met.proy$Predicted[,2],type="l",col="darkorange")
		points(proy.dates,met.proy$Predicted[,3],type="l",col="darkorange",lty=2)

	}
	# Leyenda de la figura
	par(fig=c(0,1,0,1),new=TRUE,mar=c(0,3,0,1),oma=c(0,0,0,0))
	plot(0,0,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
	legend("top",legend=c("Observado","Predicho","Proyectado","Incertidumbre en la Predicción"
			,"Incertidumbre en la Proyección")
		,pch=c(16,-1,-1,-1,-1),lty=c(-1,1,1,2,2),col=c("black","red","darkorange","red","darkorange")
		,bty="n",ncol=3,cex=0.7)
dev.off()
}

# Guardar los resultados y remover todos los archivos de la memoria
if(save.it){
	save.name	<-	paste(save.path,'Estimacion_COVID_Baq_Ciudades_',Sys.Date(),".RData",sep="")
	save(list='resultados',file=save.name)
	rm(list=ls())
}else{objects.todel	<-	ls()[!ls()%in%'resultados']
	rm(list=objects.todel)}

# Remover todos los archivos temporales excepto los resultados




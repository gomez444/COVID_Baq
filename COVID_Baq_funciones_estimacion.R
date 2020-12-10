#####################################################################################################
#
# Funciones para la estimación de parametros utilizando máxima verosimilitud via optimización numérica
# o clonación de datos.
# Ultima actualización: 9 de Diciembre de 2020
#
####################################################################################################

# Función diseñada para leer el output de getCOVID.data y modificar el formato para 
# su análisis apropiado para la estimación de parametros via Máxima Verosimilitud utilizando clonación de datos.
# x = salida de la función getCOVID.data
# ttest	= Número de observaciones para usar como datos de prueba
# K = Número de veces que se desea clonar los datos.
# location.name = nombre de la localidad para la cual se quiere dar formato a los datos. Debe ser igual a 
#				  la localidad utilizada para obtener los datos a traves de la función getCOVID.data
# Salida
# Una lista con los siguientes elementos:
# N.obs = matriz con observaciones clonadas
# lens = longitud de la serie de tiempo sin incluir la primera observación.
# K = Número de clones
# nos = vector con valores iniciales

dc.data	<-	function(x,ttest,K,location.name,data.type=c("Infect.acum","Infect.diarios","Muertos.acum","Muertos.diarios")){
	column.select<-switch(data.type,Infect.acum=1,Infect.diarios=2,Muertos.acum=3,Muertos.diarios=4)
	x			<-	x[[location.name]]
	rows.select	<-	min(which(x[,column.select]>0)):nrow(x)
	acum.infect	<- x[rows.select,column.select]
	tlength		<-	length(acum.infect)
	acum.infect.sh	<-	acum.infect[1:(tlength-ttest)]
	acum.infect.test<-	acum.infect[(tlength-ttest):tlength]
	dates	<-	as.Date(names(acum.infect))
	dates.sh<-	dates[1:(tlength-ttest)]
	dates.test	<-	dates[(tlength-ttest):tlength]

	# Definir el número de clones para la estimacion de los parametros
	K 				<- K 

	# Otros parametros necesarios 
	tslengths  <- 	tlength-ttest

	# Cloning the time series data
	N.obs=acum.infect[2:(tslengths)]
	N.obs.dc <- array(0,dim=c((tslengths-1),K)) # datos clonados sin 'no'
	for(i in 1:K){
		N.obs.dc[,i]	<-	N.obs
	}

	nos<-acum.infect[1]

	data4dc <- list(N.obs=N.obs.dc,lens = (tslengths-1), K=K, nos=nos)
	return(data4dc)
}

####################################################################################################
# Función para realizar la optimización de los parámetros del Gompertz estocástico
# usando el método de clonación de datos.
# data.list	= salida de la función dc.data
# model.type = "Simple" o "time.dependent". 
# n.chains =  número de cadenas de Markov independientes para correr.
# n.adapt = Número de iteraciones de la adaptación de las cadenas de Markov.
# n.update = Número de iteraciones para descartar de las cadenas de Markov.
# n.iter = Número de iteraciones totales en cada una de las cadenas.
# thin = Número de iteraciones a descartar entre muestras de las cadenas de Markov
# salida
# Lista con tres elementos: 
# out.dcfit = Las cadenas de Markov.
# mles.cis	= estimadores de máxima verosimilitud de los parámetros e intervalos de confianza.
# KalmanX	= Estimadores de Kalman con su intervalo de confianza para cada una de las observaciones.
dcfit <- function(data.list,model.type=c("Simple","time.dependent"),n.chains=3
					,n.adapt,n.update,n.iter,thin,inits=NULL){
	
	model	<-	switch(model.type,Simple=const.growth.model.dc,time.dependent=var.growth.model.dc)
	params	<-	switch(model.type,Simple=c("a","cc","sig")
						,time.dependent=c("beta0","beta1","cc","sig"))
	
	cl	<-	makePSOCKcluster(n.chains)
	clusterEvalQ(cl,library(dclone))
	out.jags2.dc	<-	jags.parfit(cl=cl,data.list, params=params, model=model, n.chains=n.chains
								,n.adapt=n.adapt,n.update=n.update,n.iter=n.iter,thin=thin,inits=inits)
	stopCluster(cl)

	K		   <- data.list$K
	out	 	   <- summary(out.jags2.dc)
	mles		   <- out$statistics[,1]
	mles.sds   <- out$statistics[,2]
	vars.mles <- K*(mles.sds^2)
	cis.mat	  <- cbind("2.5%"=mles-qnorm(p=0.975, mean=0,sd=1)*mles.sds
					, MLE=mles, "97.5%"=mles+qnorm(p=0.975, mean=0,sd=1)*mles.sds)
		
	sig.mle <- mles["sig"]
	cc.mles	<-	mles["cc"]
	
	switch(model.type,Simple={a.mles  <- mles["a"]}
					,time.dependent={b0.mles	<-	mles["beta0"]
									 b1.mles	<-	mles["beta1"]})

	kalman.list <- switch(model.type,Simple=list(N.obs=data.list$N.obs[,1],lens = data.list$lens
												,nos = data.list$nos,a=a.mles, cc=cc.mles
												,sigsq.inv=1/(sig.mle^2)
												)
									,time.dependent=list(N.obs=data.list$N.obs[,1],lens = data.list$lens,
														nos = data.list$nos,cc=cc.mles,
														beta1=b1.mles,beta0=b0.mles
														,sigsq.inv=1/(sig.mle^2))
														)
	kalman.params <- c("Nt")
	
	model	<-	switch(model.type,Simple=const.growth.kalman
								 ,time.dependent=var.growth.kalman)
	
	cl	<-	makePSOCKcluster(n.chains)
	clusterEvalQ(cl,library(dclone))
	kalman.jags	<-	jags.parfit(cl=cl,data=kalman.list, params=kalman.params, model=model, n.chains=n.chains
								,n.adapt=n.adapt,n.update=n.update,n.iter=n.iter,thin=thin)
	stopCluster(cl)
	out.kalman <- summary(kalman.jags)$quantiles
	KalmanX <- array(0,dim=c(data.list$len,3))
	
	KalmanX[,2]		<-	out.kalman[,"50%"]
	KalmanX[,1]		<-	out.kalman[,"2.5%"]
	KalmanX[,3]		<-	out.kalman[,"97.5%"]
	
	dimnames(KalmanX)	<-	list(rownames(data.list$N.obs),c("2.5%","MLE","97.5%"))
	
	out <- list(out.dcfit = out.jags2.dc, mles.cis = cis.mat, KalmanX=KalmanX)
	return(out)
	
}

####################################################################################################
# Función para calcular la verosimilitud negativa del modelo Gompertz deterministico
# params = vector con valores iniciales para los parametros del modelo. 
# dat = un vector con la serie de tiempo de las observaciones.
# model.type =  seleccionar el modelo para calcular la verosimilitud negativa.
# Salida
# Retorna la verosimilitud negativa del modelo gompertz deterministico
gompertz.nll	 <-function(params,dat,model.type=c("Simple","time.dependent")){
	
	switch(model.type,Simple={b0<-params[1]
					   b1<-	0
					   c <- plogis(params[2])}
			,time.dependent={b0<-params[1]
					   b1<-params[2]
					   c <-plogis(params[3])
					   })
					   	
	pred		<-	gompertz.det(b0=b0,b1=b1,c=c,len=length(dat),to=1,no=dat[1])
	loglik	<-	sum(dpois(dat[2:length(dat)],lambda=pred[2:length(dat)],log=TRUE))
	if(!is.finite(loglik)){loglik<--1*.Machine$double.xmax}
	return(-loglik)
}

####################################################################################################
# Función para minimizar la verosimilitud negativa del modelo de gompertz deterministico
# data = Serie de tiempo con el número acumulado de casos
# model.type =  Ver descripciones de los modelos. "Simple","time.dependent"
# boot = Lógico. Si se desea estimar los intervalos de confianza utilizando bootstrap paramétrico.
# nboot = número de iteraciones para calcular los intervalos de confianza con base en bootstrap paramétrico.
# inits =  Valores iniciales para la optimización.
# Salida
# Lista con:
# MLES = Estimadores de máxima verosimilitud y sus intervalos de confianza si boot = TRUE.
# loglik = log-verosimilitud del modelo.
# BIC = Criterio de información Bayesiano.

gompertz.optim<-function(data,model.type=c("Simple","time.dependent"),boot=FALSE,nboot=999,inits){
	
	if(missing(inits)){
	
		if(model.type=="Simple"){
			inits=c(log(1),qlogis(0.99))
		}else{
			inits=c(log(1),log(1),qlogis(0.99))
		}
	}
	
	optim.out	<-	optim(par=inits,fn=gompertz.nll,dat=data,model.type=model.type)
	
	mles	 <- switch(model.type,Simple=c(exp(optim.out$par[1]),plogis(optim.out$par[2]))
							,time.dependent=c(optim.out$par[1:2],plogis(optim.out$par[2]))
							)
	loglik	<-	-optim.out$value	
	BIC		<-	length(mles)*log(length(data))-2*loglik
	
	if(boot){
		mles.mat	<-	matrix(0,nrow=nboot,ncol=length(mles))
		for(i in 1:nboot){
			
			sims	<-	switch(model.type,Simple=gompertz.cond.sim(b0=mles[1],b1=0,c=mles[2],len=length(data)-1,dat=data,to=1)
											,time.dependent=gompertz.cond.sim(b0=mles[1],b1=mles[2],c=mles[3],len=length(data)-1,dat=dat,to=1)
											)
			ith.mles		<-	optim(par=inits,fn=gompertz.nll,dat=sims,model.type=model.type)
			mles.mat[i,]	<-	switch(model.type,Simple=c(exp(ith.mles$par[1]),plogis(ith.mles$par[2]))
												,time.dependent=c(ith.mles$par[1:2],plogis(ith.mles$par[2]))
												)
			
		}	
	 	mles.cis	<-	apply(mles.mat,2,quantile,prob=c(0.025,0.975))
	 	mles <- cbind(mles.cis[1,],mles,mles.cis[2,])
	 	colnames(mles)<-	c("2.5%","MLE","97.5%")
	 	rownames(mles)<-switch(model.type,Simple=c("a","c"),time.dependent=c("beta0","beta1","c"))
	}
	return(list(MLES=mles,loglik=loglik,BIC=BIC))
}

####################################################################################################
# Función para calcular la raiz cuadrada del error medio estándard del modelo logístico generalizado (Richards model).
# par = parámetros del modelo. Ver explicación en los modelos.
# x = vector con observaciones.
# Salida:
# Raiz cuadrada del error medio estándard del modelo. 

gen.log.nll	<-	function(par,x){
	gamma	<-	1/14
	K		<-	exp(par[1])
	r		<-	exp(par[2])
	theta	<-	exp(par[3])
	# beta		<-	r + gamma
	alpha	<-	plogis(par[4])
	if(alpha<0.0001){rmse<-.Machine$double.xmax
		}else{
		len	<-	length(x)-1
		predicted	<-	gen.log.mod(K,r,theta,alpha,1:len)
		rmse	<-	sqrt(mean((predicted-x[-1])^2))
	}
	return(rmse)
}

####################################################################################################
# Función para optimizar la raiz cuadrada del error medio estándard del modelo logístico generalizado (Richards model).
# x = vector con observaciones.
# inits = vector con valores iniciales de los parametros. Puede ser dejado vacío.
# hessian = Si se quiere estimar la matriz hessiana a traves de la optimización numérica.
# method = Método de optimización. Ver función optim en paquete base de R para otras opciones.
# ntest = número de observaciones para usar como prueba del modelo.
# Salida
# Una lista con:
# mles = estimadores de máxima verosimilitud y sus intervalos de confianza.
# RMSE = raiz cuadrada del error medio estándard del modelo.
# fish.inv = Matríz inversa de Información de Fisher.

gen.log.optim	<-	function(x,inits,hessian=TRUE,method="BFGS",ntest=10){
	if(missing(inits)){pars=rnorm(4)}else{pars=inits}
	opt.res	<-	tryCatch(optim(par=pars,gen.log.nll,x=x,method=method,hessian=hessian)
								,error=function(e){
								optim(par=rnorm(4),gen.log.nll,x=x,method=method,hessian=hessian)})
	rmse.vec		<-	rep(0,ntest)
	rmse.vec[1]	<-	opt.res$value
	opt.res.test	<-	list(opt.res)
	sum.val.vec	<-	0
	if(ntest>1){
		for(i in 2:ntest){
			opt.res.test[[i]]	<-	tryCatch(optim(par=rnorm(4),gen.log.nll,x=x,method=method,hessian=hessian)
											,error=function(e){
											optim(par=rnorm(4),gen.log.nll,x=x,method=method,hessian=hessian)})
			rmse.vec[i]	<-	opt.res.test[[i]]$value	
	}
		best.opt	<-	which.min(rmse.vec)
		opt.res	<-	opt.res.test[[best.opt]]
	}

	mles		<-	c(exp(opt.res$par[1:3]),plogis(opt.res$par[4]))
	names(mles)	<-	c("K","r","theta","alpha")
	out	<-	mles
	if(hessian){
		fishers.inv	<-	ginv(opt.res$hessian)
		se.mles		<-	sqrt(diag(fishers.inv))
		mles.cis		<-	c(opt.res$par-1.96*se.mles,opt.res$par+1.96*se.mles)
		mles.cis		<-	matrix(c(exp(mles.cis[1:3]),plogis(mles.cis[4]),
								exp(mles.cis[5:7]),plogis(mles.cis[8])),nrow=2
								,ncol=length(mles),byrow=TRUE)
		mles			<-	matrix(c(mles.cis[1,],mles,mles.cis[2,]),ncol=3)
		rownames(mles)<-c("K","r","theta","alpha")
		colnames(mles)<-	c("2.5%","MLE","97.5%")
		out<-list(mles=mles,RMSE=opt.res$value,fish.inv=fishers.inv)
	}
	
	return(out)
	
}

####################################################################################################
# Función para simular trayectorias del modelo logístico generalizado con base en los estimadores de 
# máxima verosimilimilitud de los parámetros.
# gen.log.optim.out = Salida de la función gen.log.optim
# len = Longitud de las trayectorias a simular.
# Bsims = Número de simulaciones.
# sim.dist = Distribución utilizada para la simulación de los datos.
# mean.width = número de observaciones antes y despues del pico para incluir en la estimación 
#				de la varianza de la distribución Negativa Binomial.
# n0		= valor inicial del número de casos nuevos
# Salida:
# Lista con:
# Predicted = Matriz con los valores predichos para el número acumulado de infectados y sus intervalos de confianza.
# Predicted.new = Matriz con los valores predichos para el número nuevo de infectados y sus intervalos de confianza.

gen.log.cis	<-	function(gen.log.optim.out,len,Bsims=1000,sim.dist=c("Poisson","NegBin"),mean.width=20,n0){
	
	mles<-gen.log.optim.out$mles[,"MLE"]
	K<-mles[1]
	r<-mles[2]
	theta<-mles[3]
	sims		<-	matrix(0,nrow=Bsims,ncol=len)
	sims.new<-	matrix(0,nrow=Bsims,ncol=len)
	v.cov	<-	diag(diag(gen.log.optim.out$fish.inv))
	mvn.draws4sim <- t(rmvnorm(n=Bsims, mean=c(log(mles[1:3]),qlogis(mles[4])),sigma=v.cov))
	mvn.draws4sim[1:3,] <-	exp(mvn.draws4sim[1:3,])
	mvn.draws4sim[4,] <-	plogis(mvn.draws4sim[4,])
	# gamma		  <-	1/14
	# beta.draws	<-	mvn.draws4sim[2,]+gamma		  
	# alpha.draws			<-	1-gamma/beta.draws
	# mvn.draws4sim	<-	rbind(mvn.draws4sim,alpha.draws)
	
	for(i in 1:Bsims){
		
		ith.K		<-	mvn.draws4sim[1,i]
		ith.r		<-	mvn.draws4sim[2,i]
		ith.theta	<-	mvn.draws4sim[3,i]
		ith.alpha	<-	min(1,mvn.draws4sim[4,i])
			
		ith.det			<-	gen.log.mod(ith.K,ith.r,ith.theta,ith.alpha,1:(len))
		if(missing(n0)){n0.hat<-ith.det[1]}else{n0.hat<-n0}
		ith.newcases.det<-	c(n0.hat,ith.det[2:len]-ith.det[1:(len-1)])
		sims.new[i,]<-	switch(sim.dist,Poisson=rpois(len,ith.newcases.det)
									 ,NegBin={	if((theta-mean.width)>length(ith.newcases.det)){
									 				mu	<-	mean(ith.newcases.det)
									 				s2	<-	var(ith.newcases.det)
									 			}else{
									 			mu		<-	mean(ith.newcases.det[(theta-mean.width):(theta+mean.width)],na.rm=TRUE)
												s2		<-	var(ith.newcases.det[(theta-mean.width):(theta+mean.width)],na.rm=TRUE)}
												if(s2<mu){s2<-mu}
												k.hat	<-	mu^2/(s2-mu)
												rnbinom(len,mu=ith.newcases.det,
															,size=k.hat)})
		if(missing(n0)){
			sims[i,]	<-	cumsum(sims.new[i,])
		}else{sims[i,]	<-	cumsum(c(ith.det[1],sims.new[i,2:len]))}

	}
	
	pred.quant		<-	t(apply(sims,2,quantile,probs=c(0.025,0.5,0.975)))
	pred.new.quant	<-	t(apply(sims.new,2,quantile,probs=c(0.025,0.5,0.975)))
	colnames(pred.quant)	<-colnames(pred.new.quant)<-	c("2.5%","50%","97.5%")
	return(list(Predicted=pred.quant,Predicted.new=pred.new.quant))
}

####################################################################################################
# Función para evaluar calcular el error medio estandar del modelo propuesto por Mellis y Littera.
# par = vector de parametros del modelo.
# x = vector de observaciones.
# Salida:
# Raiz cuadrada del error medio estándard del modelo. 

recov.nll	<-	function(par,x){
	
	c1<-exp(par[1])
	c2<-exp(par[2])
	c3<-exp(par[3])
	
	time.vec<-	1:length(x)	
	obs	<-	x
	pred<- recov.t(c1,c2,c3,time.vec)
	
	rmse	<-	sqrt(mean((pred-obs)^2))
	return(rmse)
	
}

#########################################################################################
# Estimación del Rt 
# x	= a vector con valores de nuevos infectados diarios.
# Bsims = número de simulaciones para obtener el intervalo de confianza.
# ser.int.max	=	Valor máximo del intervalo serial.
# mean.ser.int = Media de la distribución del intervalo serial.
# sd.ser.int	 = desviación estándar de la distribución del intervalo serial.
# Salida:
# Matriz con el valor promedio y los quantiles 2.5% y 97.5% de la distribución de Rt dados
# las observaciones y los parametros del intervalo serial.

rt.calc	<-	function(x,Bsims=1000,ser.int.max=20,mean.ser.int=4.79,sd.ser.int=2.9,type=c("logistic","linear")){
	
	new.cases	<-	x
	n.obs		<-	length(new.cases)
	x.vals		<-	1:n.obs
	opt			<-	switch(type,logistic=gen.log.optim(cumsum(new.cases))
							,linear=glm.nb(new.cases~x.vals))

	predicted	<-switch(type,logistic={
							tmp<-gen.log.mod(K=opt$mles["K","MLE"]
							,r=opt$mles["r","MLE"]
							,opt$mles["theta","MLE"]
							,opt$mles["alpha","MLE"],t=1:n.obs)
							c(tmp[1],tmp[2:n.obs]-tmp[1:(n.obs-1)])		
							}
							,linear=predict(opt,type="response"))

	rt.sim	<-	matrix(0,nrow=Bsims,ncol=n.obs)
	meanlog	<-	log(mean.ser.int)
	sdlog	<-	log(sd.ser.int)
		
	for(j in 1:Bsims){
		sims		<-	rpois(n.obs,predicted)
		for(i in (ser.int.max+1):n.obs){
			rt.sim[j,i]	<-	sims[i]/sum((sims[(i-ser.int.max):(i-1)]*dlnormTrunc(ser.int.max:1,meanlog,sdlog)))
		}
	}
	
	rt	<-	matrix(0,nrow=n.obs,ncol=3)
	rt[,c(1,3)]	<-	t(apply(rt.sim,2,quantile,prob=c(0.025,0.975),na.rm=TRUE))
	rt[,2]		<-	apply(rt.sim,2,mean,na.rm=TRUE)
	rownames(rt)	<-	names(new.cases)
	colnames(rt)	<-	c("lower","Mean","Upper")
	
	rt	<-	rt[-(1:(ser.int.max)),]
	return(rt)

}

###############################################################################################
# Modelos SIR (Susceptible, Infectado, Removidos)
# Función para estimar los parametros de un modelo SIR (ver modelos SIR para definición.)
# par = un vector con los valores iniciales en escala log de la tasa de contagio y el número de immigrantes
#		infectados diarios y en escala exponecial del coeficiente determinando la proporción de personas 
#		en cuarentena.
# data = una lista compuesta por los siguientes elementos ts, gamma, mu y pop. 
#		 La lista debe tener los nombres de los elementos. ts = una matriz con 3 columnas una para el número de confirmados diarios
#																una para el número de muertos y la última para el número de recuperados diarios
#															mu = tasa de mortalidad.
#															gamma = tasa de recuperación.
#															pop = población susceptible inicial.
# Salida
# log verosimilitud negativa del modelo.

sir.simp.negloglik	<-	function(par,data){
	
	ts			<-	data$ts
    n.steps		<-	ncol(ts)
	
    pars 	<- c(beta=exp(par[1]),theta=exp(par[2]),gamma =data$gamma,mu  = data$mu,b0=par[3]
    			)
	
	times	<-	1:(n.steps)
	state	<- c(S = data$pop,I = as.numeric(pars["theta"]), R = 0, C = as.numeric(pars["theta"]),D=0)
	mod.res	<-	ode(y = state, times = times, func = sir.simple, parms = pars)
		
	deaths.pred			<-	mod.res[,"D"]
	confirmed.pred		<-	mod.res[,"C"]
	recovered.pred		<-	mod.res[,"R"]
	confirmed.ll			<-	sum(dpois(ts[1,],confirmed.pred,log=TRUE))
	deaths.ll			<-	sum(dpois(ts[2,],deaths.pred,log=TRUE))
	recov.ll				<-	sum(dpois(ts[3,],recovered.pred,log=TRUE))	

	ll.pois				<-	confirmed.ll + recov.ll +deaths.ll
	
	if(!is.finite(ll.pois)){ll.pois	<-	-1*(.Machine$double.xmax)}
	
	return(-ll.pois)
		
}
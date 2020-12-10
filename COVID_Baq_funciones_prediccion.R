#####################################################################################################
#
# Funciones para la predicción de casos una vez realizadad la optimización de parametros 
# utilizando máxima verosimilitud via optimización numéricao clonación de datos.
# Ultima actualización: 9 de Diciembre de 2020
#
####################################################################################################

# Predecir las trayectorias de la epidemia con base en una conjunto de estimadores de máxima verosimilitud
# del modelo Gompertz estocástico
# dc.out = Salida de la función dcfit
# model.type = El modelo utilizado para la estimación. "Simple","time.dependent"
# tpred = Número de pasos a simular en el futuro.
# no =  valor de la primera observación que se usasra como base para la proyección al futuro.
# to = en caso que el modelo sea time.dependent, to es el número de dias despues del inicio de la epidemia
#		a partir del cual se va a comenzar con la proyección.
# init.date =  la fecha de incio de la proyección.
# Bsims = número de trayectorias a simular.
# proc.var = factor de escalamiento de la varianza del proceso.
# Salida:
# Lista con: 
# Predicted = data.frame con las fechas de la predicción y los cuantiles 25, 50 y 75 de la distribución de casos en la proyección.
# Sims = Matriz con el resultado de cada una de las trayectorias simuladas.

dc.predict	<-	function(dc.out,model.type=c("Simple","time.dependent"),tpred,no,to=1,init.date,Bsims=1000,proc.var=1){
	
	sig.mles		<-	dc.out$mles.cis["sig","MLE"]^proc.var
	if(model.type=="Simple"){
	
	a.mles	<-	dc.out$mles.cis["a","MLE"]	
	c.mles	<-	dc.out$mles.cis["cc","MLE"]
	
	}else{
		cc.mles	<-	dc.out$mles.cis["cc","MLE"]	
		b0.mles	<-	dc.out$mles.cis["beta0","MLE"]
		b1.mles	<-	dc.out$mles.cis["beta1","MLE"]
	}

	det.pred2	<-	matrix(0,nrow=Bsims,ncol=tpred)

	if(model.type=="Simple"){
		
		for(j in 1:Bsims){
		
			det.pred2[j,]	<-	const.growth.sim(parms=c(a.mles,c.mles,sig.mles),len=tpred,no=no)

		}
	
	}else{
		for(j in 1:Bsims){
			
			det.pred2[j,]	<-	var.growth.sim(parms=c(cc.mles,b0.mles,b1.mles,sig.mles),len=tpred,no=no,to=to)
		}
	}
	# If the prediction is based on the numer of new infections each day
	pred.dates	<- seq.Date(as.Date(init.date),by="days",length.out=tpred)
	pred.50		<- apply(det.pred2,2,quantile,prob=0.5,na.rm=TRUE)
	pred.25		<- apply(det.pred2,2,quantile,prob=0.25,na.rm=TRUE)	
	pred.75		<- apply(det.pred2,2,quantile,prob=0.75,na.rm=TRUE)
	colnames(det.pred2)	<-	as.character(pred.dates)
	
	out	<-	list(Predicted=data.frame(Date=pred.dates,"Q25"=pred.25,"Q50"=pred.50,"Q75"=pred.75),Sims=det.pred2)
	return(out)

}

####################################################################################################
# Función para estimar el riesgo dado el modelo Gompertz estocástico.
# dc.pred.out = salida de la función dc.predict
# threshold = valor límite para evaluar el riesgo. (Número de casos).
# Salida:
# vector con la probabilidad de sobre pasar un número de casos dado por threshold para cada día de la proyección. 
gompertz.risk	<-	function(dc.pred.out,threshold){
	
	pred.mat		<-	dc.pred.out$Sims
	thres.prob	<-	matrix(0,nrow=length(threshold),ncol=nrow(dc.pred.out$Predicted))
	colnames(thres.prob)	<-	as.character(dc.pred.out$Predicted$Date)
	rownames(thres.prob)	<-	threshold
	for(i in 1:length(threshold)){
		thres.prob[i,]	<-	apply(pred.mat,2,function(x){sum(x>=threshold[i])})/nrow(pred.mat)
	}
	return(thres.prob)	
}


mode.func	<-	function(x){
	
	 uniqx <- unique(x)
   uniqx[which.max(tabulate(match(x, uniqx)))]
	
}

####################################################################################################
# Función para predecir el resultado de un modelo Gompertz deterministico.
# optim.out = salida de la función gompertz.optim.
# model.type = "Simple", "time.dependent"
# tpred = número de días para realizar la predicción.
# no = valor inicial para la predicción.
# to = número de días a partir del inicio de la epidemia desde el que se va a hacer la proyección al futuro.
# init.date = fecha de incio de la predicción.
# Salida
# Un data frame con las fechas y los valores predichos para cada día.

gompertz.predict	<-	function(optim.out,model.type,tpred,no,to,init.date){
	
	if(class(optim.out$MLES)=="matrix"){
		params<-optim.out$MLES[,"MLE"]
	}else{
		params	<-	optim.out$MLES		
	}

	switch(model.type,Simple={b0<-log(params[1])
					   b1<-	0
					   c <- params[2]}
			,time.dependent={b0<-params[1]
					   b1<-params[2]
					   c <-params[3]
					   })
	
	pred		<-	gompertz.det(b0=b0,b1=b1,c=c,no=no,len=tpred,to=to)
	out	<-	data.frame(Date=seq.Date(init.date,by="days",length.out=tpred),Predicted=pred)
	return(out)
	
}

####################################################################################################
# Función para predecir a futuro la curva de acuerdo con el modelo logístico generalizado.
# gen.log.optim.out = salida de la función gen.log.optim
# t0 = número de días despues del inicio de la epidemia donde se quiere inciar la predicción.
# len = número de días para hacer la proyección.
# Bsims = Número de simulaciones
# sim.dist = "Poisson","NegBin". Ver descripción del modelo.
# n0		= valor inicial del número de casos nuevos
# Salida
# Lista:
# Predicted = matrix con el número acumulado de casos predichos en el futuro y los cuantules 2.5 y 97.5
# Predicted.new = matrix con el número de casos nuevos diarios predichos en el futuro y los cuantules 2.5 y 97.5
# Sims = matriz de dimensión Bsims*len con las trayectorias de las simulaciones.

gen.log.predict	<-	function(gen.log.optim.out,t0,len,Bsims=1000,sim.dist=c("Poisson","NegBin")
							,mean.width=20,n0){
	if(missing(mean.width)){mean.width=1}
	mles	<-	gen.log.optim.out$mles
	if(class(mles)=="matrix"){mles<-mles[,"MLE"]}
	K<-mles[1]
	r<-mles[2]
	theta<-mles[3]
	sims		<-	matrix(0,nrow=Bsims,ncol=len)
	sims.new<-	matrix(0,nrow=Bsims,ncol=len)
	v.cov	<-	diag(diag(gen.log.optim.out$fish.inv))
	mvn.draws4sim <- t(rmvnorm(n=Bsims, mean=c(log(mles[1:3]),qlogis(mles[4])),sigma=v.cov))
	mvn.draws4sim[1:3,] <-	exp(mvn.draws4sim[1:3,])
	mvn.draws4sim[4,] <-	plogis(mvn.draws4sim[4,])
	for(i in 1:Bsims){
		
		ith.K	<-	mvn.draws4sim[1,i]
		ith.r		<-	mvn.draws4sim[2,i]
		ith.theta	<-	mvn.draws4sim[3,i]
		ith.alpha	<-	min(1,mvn.draws4sim[4,i])
		
		ith.det	<-	gen.log.mod(ith.K,ith.r,ith.theta,ith.alpha,1:(t0+(len)))
		ith.len	<-	length(ith.det)
		if(missing(n0)){n0.hat<-ith.det[1]}else{n0.hat<-n0}
		ith.new.det	<-	c(n0.hat,ith.det[2:ith.len]-ith.det[1:(ith.len-1)])
		ith.new	<-	switch(sim.dist,Poisson=rpois(length(ith.det),ith.new.det)
									,NegBin={mu		<-	mean(ith.new.det[(theta-mean.width):(theta+mean.width)])
											s2		<-	var(ith.new.det[(theta-mean.width):(theta+mean.width)])
											if(s2<mu){s2<-mu}
											k.hat	<-	mu^2/(s2-mu)
											rnbinom(length(ith.det)
															,mu=ith.new.det
															,size=k.hat)})
		if(missing(n0)){
			sims[i,]			<-	cumsum(ith.new)[t0:(t0+len-1)]
		}else{sims[i,]<-cumsum(c(ith.det[t0],ith.new[(t0+1):(t0+len-1)]))}
		sims.new[i,]	<-	ith.new[t0:(t0+len-1)]
		
	}
	
	pred.quant	<-	t(apply(sims,2,quantile,probs=c(0.025,0.5,0.975)))
	pred.quant.new	<-	t(apply(sims.new,2,quantile,probs=c(0.025,0.5,0.975)))
	colnames(pred.quant)	<- colnames(pred.quant.new)	<-	c("2.5%","50%","97.5%")
	return(list(Predicted=pred.quant,Predicted.new=pred.quant.new,Sims=sims.new))
	
}

####################################################################################################
# Función para estimar el riesgo de sobrepasar un número determinado dado el modelo
# logistico generalizado.
# gen.log.proy.out = salida de la función gen.log.predict
# threshold = valor límite para evaluar el riesgo. (Número de casos).
# Salida
# Matriz con la probabilidad de sobrepasar threshold.
gen.log.risk	<-	function(gen.log.proy.out,threshold){
	
	pred.mat		<-	gen.log.proy.out$Sims
	thres.prob	<-	matrix(0,nrow=length(threshold),ncol=nrow(gen.log.proy.out$Predicted.new))
	rownames(thres.prob)	<-	threshold
	for(i in 1:length(threshold)){
		thres.prob[i,]	<-	apply(pred.mat,2,function(x){sum(x>=threshold[i])})/nrow(pred.mat)
	}
	return(thres.prob)	
}
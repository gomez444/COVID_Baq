# Codigo para obtener las series de tiempo para diferentes distritos en Colombia
# Distritos corresponden a los Departamentos o a un distrito de manejo especial
# wd	 = Definir el directorio para descargar o buscar los datos
# location = Ciudad, Departamento, Distrito especial o Municipio para el cual se quiere obtener los datos
# type 	   = Tipo de división administrativa para locations. "Departamento", "Municipio"
# download = lógico. Si TRUE busca los datos del INS o de la alcaldía en la web y los descarga en wd.
# date.type= Diagnostico, reporta el número acumulado y el número diario de infectados, muertos y recuperados
#			 de acuerdo con la fecha de diagnostico; Reporte, de acuerdo con la decha de Reporte en la página wb, Notificación
#			 de acuerdo con la fecha de notificación en los servicios de salud; FIS de acuerdo con la Fecha de Inicio de Sintomas.
# data.origin = Alcaldia, descarga o lee los datos que provee la alcaldia de Barranquilla. Solo funciona con type=Distritos
#				;INS descarga o lee los datos que provee el Instituto Nacional de Salud. Funciona con type= Distritos y type=Ciudad

getCOVID.data	<-	function(wd,location,type=c("Departamento","Municipio"),download,date.type=c("Diagnostico","Reporte","Notificacion","FIS")
							,data.origin=c("Alcaldia","INS")){

	switch(data.origin,Alcaldia={
			fname='casos.csv'
			# Vamos a dejar esta ruta abierta. Simplemente se debe cambiar el data.path para
			# buscar los datos provenientes de la Alcaldia con información por localidad.
			# Temporalmente hemos llamado la base de datos casos.csv pero una vez la ruta 
			# sea actualizada, la base de datos sera denominada como Alcaldia_dat_(fecha del dia).csv
			#fname		<-	paste0("Alcaldia_data_",Sys.Date(),".csv")
			#if(download){			
			#	data.path	<-	"https://www.dropbox.com/s/18gj4pw6zcyhxgn/casos.csv?dl=1"
			#	download.file(data.path,paste0(wd,"/",fname))
			#	print(paste("Data downloaded to",wd,"as",fname,"on",Sys.time()))
			#}
			#
			col.dat	<-	read.csv(paste0(wd,"/",fname),sep=";")
		},INS={
		
			fname	<-	paste0("INS_data_",Sys.Date(),".csv")

			if(download){
				# Path to Insituto Nacional de Salud Colombia and downlad data and save it with the 
				# download date in the name
				ins.path		<-	"https://www.datos.gov.co/api/views/gt2j-8ykr/rows.csv?accessType=DOWNLOAD"
				download.file(ins.path,paste0(wd,"/",fname))
				print(paste("Data downloaded to",wd,"as",fname,"on",Sys.time()))
			}
			
			col.dat	<-	read.csv(paste0(wd,"/",fname))	
		}
	)
	
	
	colnames(col.dat)	<-	iconv(colnames(col.dat),from="UTF-8",to="ASCII//TRANSLIT")
	
	if(Sys.info()[['sysname']]=="Darwin"){
		
		colnames(col.dat)	<-	gsub("'","",colnames(col.dat))
		
	}
	
	if(data.origin=="Alcaldia"){
		match.vec	<-	col.dat$LOCALIDAD
		distritos	<-	levels(match.vec)
		if(!missing(location)){
			distritos<-location
		}
	}else{
		switch(type
				,Departamento={
					match.vec	<- col.dat$Nombre.departamento	
					distritos 	<- levels(match.vec)
					if(!missing(location)){
						distritos <- location
					}
				}
				,Municipio={		
					match.vec	<-	col.dat$Nombre.municipio
					distritos 	<-	levels(col.dat$Nombre.municipio)
					if(!missing(location)){
						distritos <- location	
					}
				}
		)
	}

	switch(data.origin,Alcaldia={
								 col.dat$Fecha.de.diagnostico	<-	as.Date(col.dat$Fecha.de.diagnostico,format="%d/%m/%Y")
								 col.dat$fecha.reporte.web	<-	as.Date(col.dat$fecha.reporte.web,format="%d/%m/%Y")
								 col.dat$Fecha.de.notificacion	<-	as.Date(col.dat$Fecha.de.notificacion,format="%d/%m/%Y")
							     col.dat$Fecha.de.inicio.de.sintomas	<-	as.Date(col.dat$Fecha.de.inicio.de.sintomas,format="%d/%m/%Y")
					},INS={
							col.dat$Fecha.de.diagnostico	<-	as.Date(col.dat$Fecha.de.diagnostico,format="%d/%m/%Y")
							col.dat$fecha.reporte.web	<-	as.Date(col.dat$fecha.reporte.web,format="%d/%m/%Y")
							col.dat$Fecha.de.notificacion	<-	as.Date(col.dat$Fecha.de.notificacion,format="%d/%m/%Y")
							col.dat$Fecha.de.inicio.de.sintomas	<-	as.Date(col.dat$Fecha.de.inicio.de.sintomas,format="%d/%m/%Y")
					}
		)


	# Sumando la media de la diferencia entre la fecha de diagnostico y la fecha de notificacion a los 
	# casos que no tienen fecha de diagnostico.	
	diagnostico.nas	<-	which(is.na(col.dat$Fecha.de.diagnostico))
	mean.diag		<-	round(mean(col.dat$Fecha.de.diagnostico-col.dat$Fecha.de.notificacion,na.rm=TRUE))
	col.dat$Fecha.de.diagnostico[diagnostico.nas]	<-	col.dat$Fecha.de.notificacion[diagnostico.nas]-mean.diag

	# Sumando la media de la diferencia entre la Fecha de diagnostico y la fecha de incio de sintomas a los
	# casos asintomaticos
	FIS.nas	<-	which(is.na(col.dat$Fecha.de.inicio.de.sintomas))
	mean.FIS		<-	round(mean(col.dat$Fecha.de.diagnostico-col.dat$Fecha.de.inicio.de.sintomas,na.rm=TRUE))
	col.dat$Fecha.de.inicio.de.sintomas[FIS.nas]	<-	col.dat$Fecha.de.diagnostico[FIS.nas]-mean.FIS	
		
	column.select	<-	switch(date.type,Diagnostico=col.dat$Fecha.de.diagnostico
										,Reporte=col.dat$fecha.reporte.web
										,Notificacion=col.dat$Fecha.de.notificacion
										,FIS=col.dat$Fecha.de.inicio.de.sintomas)
	
	switch(data.origin,Alcaldia={
							col.dat$Fecha.de.muerte	<-	as.Date(col.dat$Fecha.de.muerte,format="%d/%m/%Y")
							col.dat$Fecha.de.recuperacion<-	as.Date(col.dat$Fecha.de.recuperacion,format="%d/%m/%Y")
					  },INS={
							col.dat$Fecha.de.muerte			<-	as.Date(col.dat$Fecha.de.muerte,format="%d/%m/%Y")
							col.dat$Fecha.de.recuperacion	<-	as.Date(col.dat$Fecha.de.recuperacion,format="%d/%m/%Y")		  				  	
					  }
	)
	
	time.series.list	<-	list()
	all.cases			<-	column.select
	deaths				<-	col.dat$Fecha.de.muerte[col.dat$Estado=="Fallecido"]
	recov				<-	col.dat$Fecha.de.recuperacion

	all.dates			<-	seq.Date(min(all.cases,deaths,recov,na.rm=TRUE),max(all.cases,deaths,recov,na.rm=TRUE),by="day")
	
	match.vec=iconv(match.vec,from="UTF-8",to="ASCII//TRANSLIT")
	distritos=iconv(distritos,from="UTF-8",to="ASCII//TRANSLIT")
	
	if(Sys.info()[['sysname']]=="Darwin"){
		
		match.vec	<-	gsub("'","",match.vec)
		distritos	<-	gsub("'","",distritos)
		
	}
	
	for(i in distritos){
	
		ith.distrito	<-	i
		
		if (ith.distrito!='NA')		{ith.dat			<-	col.dat[match.vec==ith.distrito,]
			}else  {ith.dat			<-	col.dat[is.na(match.vec),]}
		
		ith.column.select<-	switch(date.type,Diagnostico=ith.dat$Fecha.de.diagnostico
											,Reporte=ith.dat$fecha.reporte.web
											,Notificacion=ith.dat$Fecha.de.notificacion
											,FIS=ith.dat$Fecha.de.inicio.de.sintomas)
		casos.dates		<-	ith.column.select


		deaths.dates	<-	as.Date(ith.dat$Fecha.de.muerte[ith.dat$Estado=="Fallecido"])
		recov.dates		<-	as.Date(ith.dat$Fecha.de.recuperacion)
	
		casos.vec		<-	deaths.vec	<-	recov.vec	<-	rep(0,length(all.dates))
	
		# Infected
		casos.tmp		<-	table(casos.dates)
		casos.dates		<-	as.Date(names(casos.tmp))
		casos.vec[all.dates%in%casos.dates]	<-	casos.tmp
		cum.casos		<-	cumsum(casos.vec)

		# Deaths	
		deaths.tmp		<-	table(deaths.dates)
		if(length(deaths.tmp)==0){cum.deaths	<-	deaths.vec}else{
		deaths.dates	<-	as.Date(names(deaths.tmp))
		deaths.vec[all.dates%in%deaths.dates]		<-	deaths.tmp
		cum.deaths		<-	cumsum(deaths.vec)
		}
	
		#Recovered
		recov.tmp		<-	table(recov.dates)
		if(length(recov.tmp)==0){cum.recov	<-	recov.vec}else{
		recov.dates		<-	as.Date(names(recov.tmp))
		recov.vec[all.dates%in%recov.dates]	<-	recov.tmp
		cum.recov		<-	cumsum(recov.vec)
		}

		cum.mat		<-	matrix(c(cum.casos,casos.vec,cum.deaths,deaths.vec,cum.recov,recov.vec),nrow=length(all.dates),ncol=6
							,dimnames=list(as.character(all.dates),c("Infected_acum","Infected_daily"
																	,"Deaths_acum","Deaths_daily","Recovered_acum","Recovered_daily")))
							
		time.series.list[[i]]	<-	cum.mat
	
	}
	names(time.series.list)	<-	distritos

	return(time.series.list)

}
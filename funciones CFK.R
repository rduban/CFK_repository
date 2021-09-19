## CFA

if(!require(lavaan)){
  install.packages("lavaan")
}
require(lavaan)
mycfa<-function(mymodels,mydata){  
  #Recibe una lista que incluye los modelos y el dataframe con las variables
  #arroja una lista con los modelos especificados, objetos lavaan, la tabla de indices de ajuste,
  #y la tabla donde se determina el numero de indices de ajuste dentro los umbrales aceptables para determinar el modelo ganador
  
  if(!require(lavaan)){install.packages('lavaan')} 
  require(lavaan)       #Activar lavaan package. Se usa para Analisis de ecuaciones estructurales
  if(!require(dplyr)){install.packages('dplyr')}
  require(dplyr)        #Activar dplyr package, se utiliza para manipulacion de datos
  
  modelos<-list()
  for(i in 1:length(mymodels)){
    modelos[[i]] <-cfa(mymodels[[i]],data = mydata)  # computamos cfa y guardamos los objetos lavaan en una lista
  }
  names(modelos)<-names(mymodels)
  
  indextable0<-as.data.frame(matrix(NA,nrow=length(modelos),ncol = length(fitmeasures(modelos[[1]])))) 
  colnames(indextable0)<-names(fitmeasures(modelos[[1]]))
  row.names(indextable0)<-names(modelos)
  for(i in 1:nrow(indextable0)){
    indextable0[i,]<-round(as.vector(unlist(fitmeasures(modelos[[i]]))),3)  #obtemos los indices de ajuste de cada modelo y los guardamos en una tabla
  }
  indextable<-cbind(indextable0[3:5],indextable0[9:10],indextable0[13:14],indextable0[16],indextable0[23],indextable0[26],indextable0[29],indextable0[38:40],indextable0[42]) #seleccionamos indices de interes
  indextablevalues<-as.data.frame(matrix(NA,nrow = nrow(indextable), ncol = ncol(indextable)))
  names(indextablevalues)<-names(indextable)
  row.names(indextablevalues)<-row.names(indextable)
  indextablevalues <- rename(indextablevalues, chi2_gl = df)
  
  # a continuacion asignamos 1 cada vez que cada uno de los indices de ajuste de cada modelo se encuentran dentro de los umbrales aceptados. 
  #Se asigna 0 cuando no ocurre. Estos valores se guardan en una nueva tabla
  indextablevalues$chisq<- ifelse(indextable$chisq==min(indextable$chisq),1,0)
  indextablevalues$chi2_gl<- ifelse(nrow(mydata)>=500 & indextable$chisq/indextable$df<=2,1,ifelse(nrow(mydata)>=500 & indextable$chisq/indextable$df==min(indextable$chisq/indextable$df),1,0))
  indextablevalues$pvalue<-ifelse(indextable$pvalue<=0.05,1,0)
  indextablevalues$cfi<- ifelse(indextable$cfi>=0.9 & indextable$cfi<=1,1,0)
  indextablevalues$tli<- ifelse(indextable$tli>=0.9 & indextable$tli<=1,1,0)
  indextablevalues$nfi<- ifelse(indextable$nfi>=0.9,1,0)
  indextablevalues$pnfi<- ifelse(indextable$pnfi==max(indextable$pnfi),1,ifelse(max(indextable$pnfi)-indextable$pnfi<=0.09,1,0))
  indextablevalues$rni<-ifelse(nrow(mydata)<500 & indextable$rni==max(indextable$rni),1,0)
  indextablevalues$rmsea<- ifelse(nrow(mydata)>=100 & indextable$rmsea<=0.08,1,0)
  indextablevalues$rmsea.pvalue<-ifelse(indextable$rmsea.pvalue<=0.05,1,0)
  indextablevalues$srmr<- ifelse(nrow(mydata)<500 & indextable$srmr==min(indextable$srmr),1,0)
  indextablevalues$gfi<- ifelse(indextable$gfi==max(indextable$gfi),1,0)
  indextablevalues$agfi<- ifelse(indextable$agfi>=0.9,1,ifelse(sum(indextable$agfi>=0.9)==0 & 0.9-indextable$agfi==min(0.9-indextable$agfi),1,0))
  indextablevalues$pgfi<- ifelse(indextable$pgfi==max(indextable$pgfi),1,ifelse(sqrt((max(indextable$pgfi)-indextable$pgfi)^2)<=0.09,1,0))
  indextablevalues$ecvi<- ifelse(nrow(mydata)<500 & indextable$ecvi==min(indextable$ecvi),1,0)
  
  indextablevalues$total <- rowSums(indextablevalues) #se suman el numero de indices de ajuste dentro de los umbrales aceptados
  
  output<-NULL
  output$models<-mymodels
  output$cfamodels<-modelos
  output$all_fit_index<-indextable0
  output$fit_index<-indextable
  output$values_fit_index<-indextablevalues
  
  cat('The comparison of the models and the scoring of the fit indices were based on the guidelines of Hair, et al. (2014). Cite: Hair, J., Black, W., Babin, B. & Anderson, R. (2014). Multivariate Data Analysis (7th ed.). USA. Pearson.' )
  output
}

#calcular puntuaciones factoriales

my_scores<-function(datos,cargas){
  
  names(datos)<-paste0("q",1:ncol(datos))
  a<- datos-1
  
  b<-100*((a$q15*cargas[13,1]+a$q16*cargas[14,1]+a$q17*cargas[15,1]+a$q18*cargas[16,1]+a$q19*cargas[17,1]+a$q20*cargas[18,1]+a$q21*cargas[19,1]+a$q22*cargas[20,1]+a$q23*cargas[21,1]+a$q24*cargas[22,1])/(4*(sum(cargas[13:22,1]))))
  
  c<-100*((a$q9*cargas[8,2]+a$q10*cargas[9,2]+a$q12*cargas[10,2]+a$q13*cargas[11,2]+a$q14*cargas[12,2])/(4*(sum(cargas[8:12,2]))))
  
  d<-100*((a$q4*cargas[4,3]+a$q5*cargas[5,3]+a$q6*cargas[6,3]+a$q7*cargas[7,3])/(4*(sum(cargas[4:7,3]))))
  
  e<-100*((a$q1*cargas[1,4]+a$q2*cargas[2,4]+a$q3*cargas[3,4])/(4*(sum(cargas[1:3,4]))))
  
  out<-as.data.frame(cbind(b,c,d,e))
  apply(out,2,as.numeric)
  
  names(out)<-c("conocimientopedagC","autoeficaciaPC","comunidad","autoeficaciaTec")
  out
}

# Calificar items de conocimientos y su puntaje global

my_scores_rendim<-function(datos,grupo){
  x<-as.data.frame(matrix(NA,ncol = 1, nrow = nrow(datos)))
  names(datos)<-paste0("pc",1:ncol(datos))
  if(grupo=="inicial1"){
    
    x$pc1<- ifelse(datos$pc1=="SI" & datos$pc2=="SI",1,0)
    
    x$pc2<-ifelse(datos$pc6=="SI" & datos$pc8=="SI",1,0)
    
    x$pc3<-ifelse(is.na(datos$pc11)==TRUE,0,ifelse(datos$pc11=="3",1,0))
    
    x$pc4<-ifelse(is.na(datos$pc12)==TRUE,0,ifelse(datos$pc12=="El programa no funciona, debe capturar nuevamente el valor de la temperatura luego de encender el ventilador",1,0))
    
    x$pc5<-ifelse(is.na(datos$pc13)==TRUE,0,ifelse(datos$pc13=="6",1,0))
    
    x$pc6<-ifelse(is.na(datos$pc14)==TRUE,0,ifelse(datos$pc14=="10",1,0))
    
    x$pc7<-ifelse(is.na(datos$pc15)==TRUE,0,ifelse(datos$pc15=="II y III",1,0))
    
    ax<-strsplit(datos$pc16,split = " "); b<-c()
    for(i in 1:length(ax)){b[i]<-ifelse(ax[[i]][2]=="segunda" & ax[[i]][6]!="solamente",1,0)}
    x$pc8<-ifelse(is.na(b)==TRUE,0,b)
    
    x$pc9<-ifelse(is.na(datos$pc17)==TRUE,0,ifelse(datos$pc17=="Imagen 3",1,0))
    
    ax<-strsplit(datos$pc18,split = " "); b<-c()
    for(i in 1:length(ax)){b[i]<-ifelse(ax[[i]][6]=="120.",1,0)}
    x$pc10<-ifelse(is.na(b)==TRUE,0,b)
    
    }
  
  if(grupo=="avanzado1"){
    
    ax<-strsplit(datos$pc1,split = " "); b<-c()
    for(i in 1:length(ax)){b[i]<-ifelse(ax[[i]][1]=="alfa,",1,0)}
    x$pc1<-ifelse(is.na(b)==TRUE,0,b)
    
    ax<-strsplit(datos$pc2,split = " "); b<-c()
    for(i in 1:length(ax)){b[i]<-ifelse(ax[[i]][3]=="si",1,0)}
    x$pc2<-ifelse(is.na(b)==TRUE,0,b)
    
    x$pc3<- ifelse(is.na(datos$pc3)==TRUE,0,ifelse(datos$pc3=="Imagen 3",1,0))
    
    ax<-strsplit(datos$pc4,split = " "); b<-c()
    for(i in 1:length(ax)){b[i]<-ifelse(ax[[i]][1]=="Cada",1,0)}
    x$pc4<-ifelse(is.na(b)==TRUE,0,b)
    
    ax<-strsplit(datos$pc5,split = " "); b<-c()
    for(i in 1:length(ax)){b[i]<-ifelse(ax[[i]][6]=="r" & ax[[i]][8]!="r",1,0)}
    x$pc5<-ifelse(is.na(b)==TRUE,0,b)
    
    x$pc6<- ifelse(is.na(datos$pc6)==TRUE,0,ifelse(datos$pc6=="Cerrar el centro para carros mientras la calidad del aire sea mala, muy mala, o extremadamente mala.",1,0))
    
    ax<-strsplit(datos$pc7,split = " "); b<-c()
    for(i in 1:length(ax)){b[i]<-ifelse(ax[[i]][4]=="imprime" | ax[[i]][4]!="error",1,0)}
    x$pc7<-ifelse(is.na(b)==TRUE,0,b)
    
    x$pc8<-ifelse(is.na(datos$pc8)==TRUE,0,ifelse(datos$pc8=="Imagen 3",1,0))
  }
  
  if(grupo=="inicial2"){
    
    x$pc1<- ifelse(datos$pc1=="SI" & datos$pc2=="SI",1,0)
    
    x$pc2<-ifelse(datos$pc6=="SI" & datos$pc8=="SI",1,0)
    
    x$pc3<-ifelse(is.na(datos$pc11)==TRUE,0,ifelse(datos$pc11=="3",1,0))
    
    x$pc4<-ifelse(is.na(datos$pc12)==TRUE,0,ifelse(datos$pc12=="El programa no funciona, debe capturar nuevamente el valor de la temperatura luego de encender el ventilador",1,0))
    
    x$pc5<-ifelse(is.na(datos$pc13)==TRUE,0,ifelse(datos$pc13=="6",1,0))
    
    x$pc6<-ifelse(is.na(datos$pc14)==TRUE,0,ifelse(datos$pc14=="10",1,0))
    
    x$pc7<-ifelse(is.na(datos$pc15)==TRUE,0,ifelse(datos$pc15=="II y III",1,0))
    
    ax<-strsplit(datos$pc16,split = " "); b<-c()
    for(i in 1:length(ax)){b[i]<-ifelse(ax[[i]][2]=="segunda" & ax[[i]][6]!="solamente",1,0)}
    x$pc8<-ifelse(is.na(b)==TRUE,0,b)
    
    x$pc9<-ifelse(is.na(datos$pc17)==TRUE,0,ifelse(datos$pc17=="Imagen 3",1,0))
    
    ax<-strsplit(datos$pc18,split = " "); b<-c()
    for(i in 1:length(ax)){b[i]<-ifelse(ax[[i]][6]=="120.",1,0)}
    x$pc10<-ifelse(is.na(b)==TRUE,0,b)
    
    x$pc11<-ifelse(is.na(datos$pc19)==TRUE,0,ifelse(datos$pc19=="La botella B debe ser verde",1,0))
    
    }
  if(grupo=="avanzado2"){
    
    ax<-strsplit(datos$pc1,split = " "); b<-c()
    for(i in 1:length(ax)){b[i]<-ifelse(ax[[i]][1]=="alfa,",1,0)}
    x$pc1<-ifelse(is.na(b)==TRUE,0,b)
    
    ax<-strsplit(datos$pc2,split = " "); b<-c()
    for(i in 1:length(ax)){b[i]<-ifelse(ax[[i]][3]=="si",1,0)}
    x$pc2<-ifelse(is.na(b)==TRUE,0,b)
    
    x$pc3<- ifelse(is.na(datos$pc3)==TRUE,0,ifelse(datos$pc3=="Imagen 3",1,0))
    
    ax<-strsplit(datos$pc4,split = " "); b<-c()
    for(i in 1:length(ax)){b[i]<-ifelse(ax[[i]][1]=="Cada",1,0)}
    x$pc4<-ifelse(is.na(b)==TRUE,0,b)
    
    x$pc5<- ifelse(is.na(datos$pc5)==TRUE,0,ifelse(datos$pc5=="b",1,0))
    
    x$pc6<- ifelse(is.na(datos$pc6)==TRUE,0,ifelse(datos$pc6=="Cerrar el centro para carros mientras la calidad del aire sea mala, muy mala, o extremadamente mala.",1,0))
    
    ax<-strsplit(datos$pc7,split = " "); b<-c()
    for(i in 1:length(ax)){b[i]<-ifelse(ax[[i]][4]=="imprime" | ax[[i]][4]!="error",1,0)}
    x$pc7<-ifelse(is.na(b)==TRUE,0,b)
    
    x$pc8<-ifelse(is.na(datos$pc8)==TRUE,0,ifelse(datos$pc8=="Imagen 3",1,0))
    
    x$pc9<-ifelse(is.na(datos$pc9)==TRUE,0,ifelse(datos$pc9=="El orden de las entradas a la función no es importante. Por ejemplo, calcular(a, b, c) es lo mismo que calcular(b, a, c)",1,0))
  }
  x<-x[2:ncol(x)]
  x$total_rendim<-rowSums(x,na.rm = TRUE)
  x$total_rendim<-(x$total_rendim-mean(x$total_rendim))/sd(x$total_rendim)
  x$total_rendim<-50+x$total_rendim*10
  x
}


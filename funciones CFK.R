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
    for(i in 1:length(ax)){b[i]<-ifelse(ax[[i]][4]=="imprime" | ax[[i]][4]=="error",1,0)}
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
    for(i in 1:length(ax)){b[i]<-ifelse(ax[[i]][4]=="imprime" | ax[[i]][4]=="error",1,0)}
    x$pc7<-ifelse(is.na(b)==TRUE,0,b)
    
    x$pc8<-ifelse(is.na(datos$pc8)==TRUE,0,ifelse(datos$pc8=="Imagen 3",1,0))
   
    ax<-strsplit(datos$pc9,split = " "); b<-c()
    for(i in 1:length(ax)){b[i]<-ifelse(ax[[i]][1]=="x::3,",1,0)}
    x$pc9<-ifelse(is.na(b)==TRUE,0,b)
  }
  x<-x[2:ncol(x)]
  x$total_rendim<-rowSums(x,na.rm = TRUE)
  x$total_rendim<-(x$total_rendim-mean(x$total_rendim))/sd(x$total_rendim)
  x$total_rendim<-50+x$total_rendim*10
  x
}

my_t.test<-function(varscon,varfact,paired,var.equal){
  varfact<-factor(varfact)
  tabla<-as.data.frame(matrix(NA,ncol = 10, nrow = 1)); names(tabla)<-c("t value","df","p.value","d.cohen",paste("M",levels(varfact)[1]),paste("SD",levels(varfact)[1]),paste("n",levels(varfact)[1]),paste("M",levels(varfact)[2]),paste("SD",levels(varfact)[2]),paste("n",levels(varfact)[2]))
  a<-tapply(varscon,varfact,FUN = function(x){round(na.omit(mean(x)),2)})
  b<-tapply(varscon,varfact,FUN = function(x){round(na.omit(sd(x)),2)})
  n<-tapply(varscon,varfact,FUN = function(x){round(na.omit(length(x)),2)})
  
  if(paired==TRUE){
      x<-t.test(varscon~varfact, paired = TRUE, var.equal = ifelse(var.equal==TRUE,TRUE,FALSE))
      tabla[1,4]<-sqrt((round(x$estimate/min(tapply(varscon,varfact,sd)),3))^2)
  }
  if(paired==FALSE){
    x<-t.test(varscon~varfact, paired = FALSE,var.equal = ifelse(var.equal==TRUE,TRUE,FALSE))
    tabla[1,4]<-sqrt((round((max(x$estimate)-min(x$estimate))/mean(tapply(varscon,varfact,sd)),3))^2)
  }
  
  tabla[1,1]<- round(x$statistic,2)
  tabla[1,2]<- x$parameter
  tabla[1,3]<- round(x$p.value,3)
  tabla[1,5]<-a[1]
  tabla[1,6]<-b[1]
  tabla[1,7]<-n[1]
  tabla[1,8]<-a[2]
  tabla[1,9]<-b[2]
  tabla[1,10]<-n[2]
  
  tabla
  }

my_anova2x2<-function(modelos,level){
#a<-aov(lm(autoeficaciaTec ~fase*genero+fase*contexto+fase*niv.ensenanza+fase*area, data = b))
tabla<-as.data.frame(matrix(NA,ncol = 5, nrow = length(modelos)))
dimnames(tabla)<-list(names(modelos),c("df","df.denom","F value","p value","eta.sq"))
for(i in 1:length(modelos)){
  a<-summary(modelos[[i]])
  z<-etaSquared(modelos[[i]],type=2,anova = FALSE)
  
  if(level==3){
   b<-a[[1]][3,]  
   df_dn<-a[[1]][4,1]
   eta.sq<-z[3,1]  
  }
  if(level==2){
   b<-a[[1]][2,]  
   df_dn<-a[[1]][3,1]
   eta.sq<-z[2,1]
  }
  if(level==1){
   b<-a[[1]][1,]  
   df_dn<-a[[1]][2,1]
   eta.sq<-z[1,1]
  }
  tabla[i,1]<-b[1]; tabla[i,2]<-df_dn; tabla[i,3:4]<-round(b[4:5],3);tabla[i,5]<-round(eta.sq,3)
  }
tabla
}

my_aciertos<-function(datos ,var_grupo, titulo, curso){
  # datos = dataframe
  # var_grupo = vector con valores dicotomicos
  # nombre del vector dicotomico para que aparezca en el titulo de la grafica
  if(is.factor({{var_grupo}})==TRUE){
    x<-as.vector(levels({{var_grupo}})); x<-x[!is.na(x)]
  }else{x<-as.vector(unique({{var_grupo}})); x<-x[!is.na(x)]}
  
  etiquetas<-paste0("desempenio",1:length(x))
  lista<-NULL
  if(curso=="inicial1"){
    aciertos<-as.data.frame(matrix(NA,ncol = length(x), nrow =10))
    a<-as.data.frame(matrix(NA,ncol = length(x), nrow = 10))
    for(i in 1:length(x)){
      lista[[i]]<- datos[9:18] %>% filter({{var_grupo}}==x[i])
      j<- round(100*colSums(lista[[i]])/nrow(lista[[i]]),2)
      aciertos[i]<-j
      a[i]<-rep(x[i],10)
    }
    Categorias<-rep(c("Problemas de caso","Problemas de caso","Diagramas de flujo","Diagramas de flujo","Operaciones internas","Condicionales","Lógica booleana","Lógica booleana","Condicionales","Bucles"),length(x))
    niveles<-c(paste0("pc",1:10))
  }
  if(curso=="avanzado1"){
    aciertos<-as.data.frame(matrix(NA,ncol = length(x), nrow =8))
    a<-as.data.frame(matrix(NA,ncol = length(x), nrow = 8))
    for(i in 1:length(x)){
      lista[[i]]<- datos[9:16] %>% filter({{var_grupo}}==x[i])
      j<- round(100*colSums(lista[[i]])/nrow(lista[[i]]),2)
      aciertos[i]<-j
      a[i]<-rep(x[i],8)
    }
    Categorias<-rep(c("Operaciones internas","Condicionales pseudo-código","Condicionales MakeCode","Condicionales MakeCode","Operaciones internas","Bucles","Condicionales pseudo-código","Condicionales MakeCode"),length(x))
    niveles<-c(paste0("pc",1:8))
  }
  if(curso=="inicial2"){
    aciertos<-as.data.frame(matrix(NA,ncol = length(x), nrow =11))
    a<-as.data.frame(matrix(NA,ncol = length(x), nrow = 11))
    for(i in 1:length(x)){
      lista[[i]]<- datos[9:19] %>% filter({{var_grupo}}==x[i])
      j<- round(100*colSums(lista[[i]])/nrow(lista[[i]]),2)
      aciertos[i]<-j
      a[i]<-rep(x[i],11)
    }
    Categorias<-rep(c("Problemas de caso","Problemas de caso","Diagramas de flujo","Diagramas de flujo","Operaciones internas","Condicionales","Lógica booleana","Lógica booleana","Condicionales","Bucles","Lógica booleana"),length(x))
    niveles<-c(paste0("pc",1:11))
  }
  if(curso=="avanzado2"){
    aciertos<-as.data.frame(matrix(NA,ncol = length(x), nrow =9))
    a<-as.data.frame(matrix(NA,ncol = length(x), nrow = 9))
    for(i in 1:length(x)){
      lista[[i]]<- datos[9:17] %>% filter({{var_grupo}}==x[i])
      j<- round(100*colSums(lista[[i]])/nrow(lista[[i]]),2)
      aciertos[i]<-j
      a[i]<-rep(x[i],9)
    }
    Categorias<-rep(c("Operaciones internas","Condicionales","Condicionales","Condicionales","Operaciones internas","Bucles","Funciones","Condicionales","Funciones"),length(x))
    niveles<-c(paste0("pc",1:9))
  }
  desempenio<-f.reubica(aciertos)
  preguntas<-rep(paste0("pc",1:nrow(aciertos)),length(x))
  grupos<-f.reubica(a)
  aciertos<- as.data.frame(cbind(desempenio,preguntas,Categorias,grupos))
  colnames(aciertos)<- c("Porcentaje", "Pregunta","Categoria","Fase")
  aciertos$Porcentaje<- as.numeric(as.character(aciertos$Porcentaje))
  
  aciertos$Pregunta<-factor(aciertos$Pregunta, levels = niveles)
  aciertos$Categoria<- factor(aciertos$Categoria, levels= unique(Categorias))
  aciertos$Fase<-factor(aciertos$Fase, levels = x)
  
  ggplot(data=aciertos, aes(x=Pregunta, y=Porcentaje, fill=Fase)) +
    geom_bar(stat="identity", position = "dodge")+ scale_fill_manual(values=brewer.pal(n = length(unique(Categorias)), name = "Pastel2"))+ geom_text(aes(label=paste(Porcentaje, "%", sep="")), size = 2.7, vjust=0, color="black", position = position_dodge(0.9), size=3.1)+(labs(title = paste("Aciertos por",titulo), x = "", y = "Porcentaje de aciertos"))+
    ylim(0,100) + facet_wrap(aciertos$Categoria, scales = "free")
}

my_aciertos_gr<-function(datos ,var_grupo, titulo, curso){
  # datos = dataframe
  # var_grupo = vector con valores dicotomicos
  # nombre del vector dicotomico para que aparezca en el titulo de la grafica
  if(is.factor(var_grupo)==TRUE){
    x<-as.vector(levels(var_grupo)); x<-x[!is.na(x)]
  }else{x<-as.vector(unique(var_grupo)); x<-x[!is.na(x)]}
  
  etiquetas<-paste0("desempenio",1:length(x))
  lista<-NULL
  if(curso=="inicial1"){
    aciertos<-as.data.frame(matrix(NA,ncol = length(x), nrow =10))
    a<-as.data.frame(matrix(NA,ncol = length(x), nrow = 10))
    for(i in 1:length(x)){
      lista[[i]]<- datos[9:18] %>% filter(var_grupo==x[i])
      j<- round(100*colSums(lista[[i]])/nrow(lista[[i]]),2)
      aciertos[i]<-j
      a[i]<-rep(x[i],10)
    }
    Categorias<-rep(c("Problemas de caso","Problemas de caso","Diagramas de flujo","Diagramas de flujo","Operaciones internas","Condicionales","Logica booleana","Logica booleana","Condicionales","Bucles"),length(x))
    niveles<-c(paste0("pc",1:10))
  }
  if(curso=="avanzado1"){
    aciertos<-as.data.frame(matrix(NA,ncol = length(x), nrow =8))
    a<-as.data.frame(matrix(NA,ncol = length(x), nrow = 8))
    for(i in 1:length(x)){
      lista[[i]]<- datos[9:16] %>% filter(var_grupo==x[i])
      j<- round(100*colSums(lista[[i]])/nrow(lista[[i]]),2)
      aciertos[i]<-j
      a[i]<-rep(x[i],8)
    }
    Categorias<-rep(c("Operaciones internas","Condicionales pseudo-codigo","Condicionales MakeCode","Condicionales MakeCode","Operaciones internas","Bucles","Condicionales pseudo-codigo","Condicionales MakeCode"),length(x))
    niveles<-c(paste0("pc",1:8))
  }
  if(curso=="inicial2"){
    aciertos<-as.data.frame(matrix(NA,ncol = length(x), nrow =11))
    a<-as.data.frame(matrix(NA,ncol = length(x), nrow = 11))
    for(i in 1:length(x)){
      lista[[i]]<- datos[9:19] %>% filter(var_grupo==x[i])
      j<- round(100*colSums(lista[[i]])/nrow(lista[[i]]),2)
      aciertos[i]<-j
      a[i]<-rep(x[i],11)
    }
    Categorias<-rep(c("Problemas de caso","Problemas de caso","Diagramas de flujo","Diagramas de flujo","Operaciones internas","Condicionales","Logica booleana","Logica booleana","Condicionales","Bucles","Lógica booleana"),length(x))
    niveles<-c(paste0("pc",1:11))
  }
  if(curso=="avanzado2"){
    aciertos<-as.data.frame(matrix(NA,ncol = length(x), nrow =9))
    a<-as.data.frame(matrix(NA,ncol = length(x), nrow = 9))
    for(i in 1:length(x)){
      lista[[i]]<- datos[9:17] %>% filter(var_grupo==x[i])
      j<- round(100*colSums(lista[[i]])/nrow(lista[[i]]),2)
      aciertos[i]<-j
      a[i]<-rep(x[i],9)
    }
    Categorias<-rep(c("Operaciones internas","Condicionales","Condicionales","Condicionales","Operaciones internas","Bucles","Funciones","Condicionales","Funciones"),length(x))
    niveles<-c(paste0("pc",1:9))
  }
  desempenio<-f.reubica(aciertos)
  preguntas<-rep(paste0("pc",1:nrow(aciertos)),length(x))
  Fase<-f.reubica(a)
  aciertos<- as.data.frame(cbind(desempenio,preguntas,Categorias,Fase))
  colnames(aciertos)<- c("Porcentaje", "Pregunta","Categoria","Fase")
  aciertos$Porcentaje<- as.numeric(as.character(aciertos$Porcentaje))
  
  aciertos$Pregunta<-factor(aciertos$Pregunta, levels = niveles)
  aciertos$Categoria<- factor(aciertos$Categoria, levels= unique(Categorias))
  aciertos$Fase<-factor(aciertos$Fase, levels = c("pre-test","post-test"))
  
  ggplot(data=aciertos, aes(x=Pregunta, y=Porcentaje, fill=Fase)) +
    geom_bar(stat="identity", position = "dodge")+ scale_fill_manual(values=brewer.pal(n = length(unique(Categorias)), name = "Pastel2"))+ geom_text(aes(label=paste(Porcentaje, "%", sep="")), size = 2.7, vjust=0, color="black", position = position_dodge(0.9), size=3.1)+(labs(title = paste("Aciertos de",titulo), x = "", y = "Porcentaje de aciertos"))+
    ylim(0,100) + facet_wrap(aciertos$Categoria, scales = "free")
}

my_graph.practicas<-function(datos,grupo, fase){
theme_set(theme_pubclean())
datos<-na.omit(datos)
diamonds.frac<-dplyr::sample_frac(datos)

if(grupo==1){
names(diamonds.frac)<-c("Actividades desconectads","Usa-modifica-crea","Clase magistral","Enseñanza explicita","Marcha silenciosa", "Apren. basado en proyectos")}
if(grupo==2){
names(diamonds.frac)<-c("Dar respuesta correcta","Sugerir ir paso a paso","Sugerir revisar notas","Sugerir revisar memorias","Sugerir volver a leer","Usar diferentes valores","Volver a explicar")}
######
variables <- colnames(diamonds.frac)
diamonds.frac[ diamonds.frac == 0 ] <- 'No lo conozco'

diamonds.m <- reshape::melt(diamonds.frac,id.vars=0)

diamonds.m$value<-factor(diamonds.m$value,levels = c("No lo conozco","1","2","3","4","5","6","7","8","9","10"))

ggplot(diamonds.m ,aes(x=variable, y=value, color = value))+
  geom_jitter(size=1)+
  labs(title="Preferencia practica pedagogica", subtitle = fase,x= "Practica pedagogica",y = "Nivel de preferencia")+
  theme(axis.text.x = element_text(face = "italic", angle = 0, hjust = 0.5, vjust = 0.5, size=rel(1)))+
   scale_color_brewer(palette="Set3")+theme(legend.position="none")
  
  }

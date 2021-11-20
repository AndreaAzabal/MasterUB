library(rgdal) 
library(rgeos)
library(data.table)
library(osmdata)
library(maptools)
library(leaflet)
library(sp)
library(dplyr)
library(geosphere)
library(pander)
library(ggcorrplot)
library(ggplot2)
library(lmtest) 
library(fBasics) 
library(MASS)
library(earth)
library(fitdistrplus)
library(mfx)
library(pROC)
library(glmnet)
library(spdep)
library(spgwr)
library(astsa)
library(forecast)
library(car)
library(pls)
library(MPV) 

#Funcion para establecer el percentil de una determinada variable
perc.rank <- function(x) trunc(rank(x))/length(x)

#Función para descargar puntos de OSM
Descarga_OSM<-function(ciudad="Madrid, Spain",key='building',value = "hospital"){

#Descargo la Iformación
mapa1 <- opq(bbox = ciudad)
Poligonos_dentro <- add_osm_feature(mapa1, key = key, value = value)
df <- osmdata_sp(Poligonos_dentro)
#Centroides de cada polígono + representación
spChFIDs(df$osm_polygons) <- 1:nrow(df$osm_polygons@data)
centroides <- gCentroid(df$osm_polygons, byid = TRUE)
names<-df$osm_polygons$name

#Creo los Buffers de Hospitales. En menos de 200 metros. 
buffer <- gBuffer(centroides, byid = TRUE, width = 0.002)

#Convierto en Spatial Polygon DataFrame
buffer <- SpatialPolygonsDataFrame(buffer, data.frame(row.names = names(buffer), n = 1:length(buffer)))
#Combino los Polígonos que se entrecruzan
gt <- gIntersects(buffer, byid = TRUE, returnDense = FALSE)
ut <- unique(gt); nth <- 1:length(ut); buffer$n <- 1:nrow(buffer); buffer$nth <- NA
for(i in 1:length(ut)){
  x <- ut[[i]];  buffer$nth[x] <- i}
buffdis <- gUnaryUnion(buffer, buffer$nth)

#Combino los Polígonos que se entrecruzan otra vez.
gt <- gIntersects(buffdis, byid = TRUE, returnDense = FALSE)
ut <- unique(gt); nth <- 1:length(ut)
buffdis <- SpatialPolygonsDataFrame(buffdis, data.frame(row.names = names(buffdis), n = 1:length(buffdis)))
buffdis$nth <- NA
for(i in 1:length(ut)){
  x <- ut[[i]];  buffdis$nth[x] <- i}
buffdis <- gUnaryUnion(buffdis, buffdis$nth)

sd<-list(centroides,buffdis,names)

return(sd)

}

#Crea histograma y representa como se comporta la variable respuesta en el
Hist<-function(bbdd_fff,response,predicted,var,n,breaks=8){
  
  names<-colnames(bbdd_fff)
  bbdd_fff$predicted<-predicted
  bbdd_fff$response<-response
  bbdd_fff$VIVO<-1
  q1<- bbdd_fff %>% 
    group_by(cut(var,breaks=breaks)) %>% 
    summarise(Exposicion = sum(VIVO),
              Frecuencia = sum(response)/sum(VIVO),
              Predicted = sum(predicted)/sum(VIVO)
              ) 
 
  q1<-as.data.table(q1)
  
  c1<-q1[,1]
  q1$c1<-q1[,1]
  
  ff<-ggplot(q1,aes(x=c1)) + 
    geom_bar(aes(y=(Exposicion/sum(Exposicion))*mean(Frecuencia)*2, fill="% de Exposicion"), stat = "identity")+ 
    geom_point(aes(y=Frecuencia, colour="Frecuencia"), group=1)+
    geom_line(aes(y=Frecuencia, colour="Frecuencia"), group=1)+
    geom_point(aes(y=Predicted, colour="Predicted"), group=1)+
    geom_line(aes(y=Predicted, colour="Predicted"), group=1)+    
    xlab("var") + ylab("") +
    #scale_y_continuous(limits = c(0, 0.45),breaks = c(0:45/100))+
    theme(legend.key=element_blank(),legend.title=element_blank(),
          legend.box="horizontal", legend.position = "bottom",
          axis.text.x = element_text(angle = 60)) +
    ggtitle(names[n])
  
  return(ff)
  
  
}

#Crea diagrama de cajas
diag_cajas<-function(table,filas=6,columnas=3,horizontal=TRUE){
  
  n<-ncol(table)
  mar.default <- c(2,2,2,2) 
  #par(mar=c(1,1,1,1))
  par(mfrow = c(filas, columnas), mar = mar.default + c(2, 0, 0, 0)) 

  sapply(seq(1,length(table)),function(j)boxplot(table[,j],main=colnames(table)[j],horizontal = horizontal,xlab="",col="blue"))

  
  
  
}

#Pinta a Nivel Punto una base de datos y Variable COlor y Tama?o
pl_pt<-function(df,size2,color2,dd=5,sz=500){
  
  volterars=0
  volterarc=0
  
  if (!is.numeric(size2)) {  df$size<-as.numeric(as.factor(size2)) }
  if (!is.numeric(color2)) { df$color<-as.numeric(as.factor(color2))}
  if (is.numeric(size2)) {  df$size<-(size2) }
  if (is.numeric(color2)) { df$color<-(color2)}
  x<-dd 
  dd<-seq(0,1,1/dd)
  
  if (volterars==1){      df$size<-(max(df$size)+1-df$size)    }
  if (volterarc==1){      df$color<-(max(df$color)+1-df$color)    } 
  
  
  if (length(unique(df$color))<10){    pal <- colorBin(palette = "RdYlBu", domain = df$color ,bins = length(levels(as.factor(df$color))) , na.color = "grey40", reverse = T) }
  if (length(unique(df$color))>=10){   pal <- colorBin(palette = "RdYlBu", domain = df$color ,bins = unique(quantile(df$color, dd )), na.color = "grey40", reverse = T) }
  
  a<-as.character(cut(as.numeric(as.factor(df$size)),breaks=x))
  a<-as.numeric(as.factor(a))
  
  
   pintar<-leaflet() %>%
    addTiles() %>%
    addLegend(pal = pal, values = round(df$color, 1), position = "bottomright", title = "") %>%
    addCircles(data=df,lng =df$LONG ,lat =df$LAT , stroke = FALSE, opacity = 0.5,fillOpacity = 0.5,
               color =pal(df$color),radius=a*sz)
  
  return(pintar)
  
  
  
}



## -------------------------------------------------------------------------
# nombres de los paquetes de Bioconductor que vamos a instalar
bioc.p=c("DESeq2", "vsn", "apeglm", "genefilter", "IHW", "edgeR")


## ---- eval=F--------------------------------------------------------------
## 
## # primero se instala el paquete BiocManager que usaremos para descargar otros paquetes
## # (sí, R es un caos de paquetes)
install.packages("BiocManager")
## 
## # descargamos los paquetes que espcificamos anteriormente en la lista bioc.p
BiocManager::install(version = "3.12", force = T)
BiocManager::install(bioc.p, force = T)
## 




## ---- eval=T--------------------------------------------------------------

cran.p=c("dplyr", # manipulación de datos
         "tidyr", # manipulación de datos
         "ggplot2", # generación de gráficos
         "pheatmap", # gráficos de mapas de calor
         "RColorBrewer", # paletas de colores para los gráficos
         "PoiClaClu", # cálculo de distancias de Poisson
         #"glmpca",
         "ggbeeswarm", 
         'gridExtra', 
         'colorspace')



## ---- eval=F--------------------------------------------------------------
## 
new.packages=cran.p[!(cran.p %in% installed.packages()[,"Package"])]

if (length(new.packages)>0) {
  install.packages(new.packages)
}



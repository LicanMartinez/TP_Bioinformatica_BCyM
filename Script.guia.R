## ----Carga de librerías, warning=FALSE, message=FALSE, results='hide'-----

# Definimos la lista de nombres de los paquetes que R va a llamar

pckg.ls=c("DESeq2", "vsn", "apeglm", "genefilter", "IHW", "edgeR", # análiis de datos moleculares (bioconductor)
         "dplyr", "tidyr", # manipulación de datos
         "ggplot2", # generación de gráficos
         "pheatmap", # gráficos de mapas de calor
         "RColorBrewer", # paletas de colores para los gráficos
         # "PoiClaClu", # cálculo de distancias de Poisson
         # "ggbeeswarm", 
         'gridExtra', # múltiplies gráficos en un único panel
         'colorspace' # configuraciones de color
         )

# E importamos los paquetes de esa lista

lapply(pckg.ls, # para cada nombre de la lista de paquetes instalados...
       require, # aplicar la función "require" para importarlo
       quietly =T,
       character.only = TRUE) # (ignorar este argumento, es para las mañas de R)



## -------------------------------------------------------------------------
# importación
cts <- read.delim("./Conteos.tsv", row.names = "gene")



## -------------------------------------------------------------------------
# filas de 1 a 6, columnas de 1 a 3
cts[1:6, 1:3]





## -------------------------------------------------------------------------

samples <- read.csv("./Muestras_reord.tsv", row.names=1)



## -------------------------------------------------------------------------
# filas de 1 a 8, todas las columnas
samples[1:8,]



## -------------------------------------------------------------------------

# nombres de las categorías a conservar (todas en este caso)
tejidos=c('CHE', 'THE', 'BRN', 'GEN')
sex=c('F', 'M')
trat=c('ctr', 'dsx')


# Índice T/F de muestras que pertenecen a las combinaciones de categorías especificadas
#   El operador & indica T sólo si se cumplen las condiciones a izquierda y derecha a la vez 
#   (T&T = T; T&F = F; F&F = F)
idx.muestras= samples$Tissue %in% tejidos & # TRUE para las muestras pertenecientes (%in%) a los tejidos defeinidos (todos)
  samples$Sex %in% sex & # TRUE para las muestras pertenecientes (%in%) a los sexos defeinidos (todos)
  samples$Treatment %in% trat # TRUE para las muestras pertenecientes (%in%) a los tratamientos defeinidos (todas)



## ---- echo=T, results='hide'----------------------------------------------

cbind(samples, idx.muestras) 



## -------------------------------------------------------------------------

cts=cts[,rownames(samples)] # para corregir diferencias en el orden de las columnas-filas entre tablas

all(rownames(samples) == colnames(cts)) # debe devolver TRUE



## ---- warning=F-----------------------------------------------------------
#objeto DESeqDataSet
data_0 <- DESeqDataSetFromMatrix(countData = cts[,idx.muestras],
                                 colData = samples[idx.muestras, c("Sex", "Tissue", "Treatment")],
                                 design = ~  Sex + Tissue + Treatment) 


## -------------------------------------------------------------------------
data_0


## -------------------------------------------------------------------------
# Probar con diferentes valores
reads.minimos=100

# genes sobre el umbral (suficientes reads)
c(rowSums(counts(data_0)) >= reads.minimos ) %>% sum()



## -------------------------------------------------------------------------
reads.minimos=10

data_filtered <- data_0[rowSums(counts(data_0)) >= reads.minimos,]

# Cantidad de genes finales
n.genes=nrow(data_filtered)



## -------------------------------------------------------------------------
data_var.st=vst(data_filtered, blind = FALSE)


## -------------------------------------------------------------------------

# cálculo de las nuevas dimensiones (proyección sobre dos componentes principales)
pcaData <- plotPCA(data_var.st, 
                   intgroup = c("Tissue", "Sex", "Treatment"), 
                   returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData$tt=paste0(pcaData$Tissue, '-', pcaData$Treatment)



## -------------------------------------------------------------------------
# Configuración de los colores de los puntos 
colores=c(BRN='dodgerblue', 
               CHE='limegreen', 
               GEN='darkgoldenrod2', 
               THE='coral1')[tejidos]

colores.trat=colorspace::lighten(colores, 0.5)
names(colores.trat)=paste0(tejidos, '-dsx')

colores.ctrl=colorspace::darken(colores.trat, 0.4)
names(colores.ctrl)=paste0(tejidos, '-ctr')

colores2=c(colores.trat, colores.ctrl)


## ---- fig.height=5, fig.width=6, dev='png', dpi=300, fig.align='center', results='hide'----

# Generación inicial del gráfico
ggplot(pcaData, # datos de entrada
       
       # aspectos del gráfico a ser mapeados a valores de variables presentes en los datos
       aes(x = PC1, # posición en X
           y = PC2, # posición en Y
           fill = tt)) + # color de relleno
  
  # Recién acá le indicamos que queremos dibujar puntos sobre el gráfico.
  geom_point(size =2.5, shape=23) + 
  # Sus coordenadas y colores van a ser definidos por el mapeo especificado en aes().
  # Tamaño y forma son especificados fuera de aes() ya que, en este caso, 
  #     no dependen de variables de los datos
  
  # Especificación del conjunto de colores y su leyenda
  scale_fill_manual(values = colores2, 
                    name=' ctr    dsx', 
                    labels=c(rep('', length(tejidos)), tejidos))+
  
  guides(fill=guide_legend(ncol=2)) +
  
  # titulos de ejes
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  
  facet_grid(Sex ~ .)+ # ¿?
  
  theme_bw() # estética



## -------------------------------------------------------------------------
# cálculo de distancias
sampleDists <- dist(t(assay(data_var.st)))

# transformación del tipo de objeto R (class: dist -> matrix)
sampleDistMatrix <- as.matrix( sampleDists )

# homologar nombres de filas y columnas
rownames(sampleDistMatrix) <- rownames(samples)[idx.muestras]
colnames(sampleDistMatrix) <- rownames(samples)[idx.muestras]

sampleDistMatrix[1:4, 1:4] # ver primeras 4 filas y 4 columnas


## ---- fig.height=6.8, fig.width=8, dev='png', dpi=400, fig.align='center', results='hide'----
# çonfiguración de la paleta de colores
heatmap.1.colors <- colorRampPalette( brewer.pal(9, "YlOrRd")[6:1] )(255)

pheatmap(sampleDistMatrix,
         col = heatmap.1.colors, border_color = NA, 
         fontsize_row = 4.5, fontsize = 11, # tamaño de las letras
         cluster_cols = F, cluster_rows = F, # mantener el orden de filas y columnas
         annotation_col  = samples[idx.muestras,1:3], # información para el código de anotación por colores de columnas
         annotation_row  = samples[idx.muestras,1:3], # información para el código de anotación por colores de filas 
         show_colnames=F, show_rownames=T, # nombres de muestras en filas
         annotation_names_row=F, legend = F, # mostrar qué significa cada color
         
         # configuración de los colores de cada grupo de muestras
         annotation_colors = list(Sex=c(M='slateblue2', 
                                        F='hotpink1')[sex], 
                                  
                                  Treatment=c(ctr='gray90', 
                                              dsx='gray10')[trat], 
                                  
                                  Tissue=c(BRN='dodgerblue', 
                                           CHE='limegreen', 
                                           GEN='darkgoldenrod2', 
                                           THE='coral1')[tejidos]
                                  )
         )



## ---- warning=F-----------------------------------------------------------

# nombres de las categorías a conservar
tejidos=c('CHE') # tejido cefálico
sex=c('F', 'M') # ambos sexos
trat=c('ctr') # individuos control

# índice T/F de muestras que pertenecen a las combinaciones de categorías especificadas
# el operador & indica T sólo si se cumplen las condiciones a izquierda y derecha a la vez 
# (T&T = T; T&F = F; F&F = F)
idx.muestras= samples$Tissue %in% tejidos &
  samples$Sex %in% sex &
  samples$Treatment %in% trat

# objeto DESeqDataSet (conteos + información de muestras)
data.DE_0 <- DESeqDataSetFromMatrix(countData = cts[,idx.muestras],
                                       colData = samples[idx.muestras, 
                                                         c("Sex", "Tissue", "Treatment")],
                                       design = ~  Sex)


reads.minimos=10

data_DE <- data.DE_0[rowSums(counts(data.DE_0)) >= reads.minimos,]

# Cantidad de genes incluidos en el análisis
nrow(data_DE)



## ---- results='hide'------------------------------------------------------
de.an=DESeq(data_DE)


## -------------------------------------------------------------------------
res = results(de.an, contrast = c('Sex', 'M', 'F')) %>% 
  as.data.frame() %>% 
  dplyr::select(-c(4:5)) # borrar algunas columans que no son de interés

head(res, 5)



## ---- results='hide'------------------------------------------------------

res2=res[order(abs(res$log2FoldChange), decreasing = T),] %>% # re-ordenamiento
  as.data.frame() %>% # llevar a formato (class) data.frame
  filter(padj<=0.05) # filtrar por significancia (p-valor < 0.05)

# print primeras filas
print(head(res2))



## -------------------------------------------------------------------------

write.table(res2%>% # tabla completa
              mutate(gene=rownames(res2)) %>% # agregamos columna de genes
              mutate(genomic.sequence=NA) %>% # columna vacía para llenar después
              head(), 
            file = './head.results_samples.CHE.ctrl_DE.MvsF.csv', 
            sep = '\t', row.names = F)



## -------------------------------------------------------------------------

# se pasa a clase "matrix" y se transpone la tabla (se acuesta) para normalizar por columnas
cts.trans=cts %>% 
  as.matrix() %>% 
  t() 

# vector de tamaños de librerías (en unidades de 10 millones)
norm.vec=(samples$lib.size/1e7)

# normalización por columnas (dividir la matriz por el vector)
# Nota: en R esto es al reves que en las operaciones matriciales de algebra, 
#       donde los vectores se aplican sobre las filas
cts.trans.norm=cts.trans/norm.vec

# vuelvo a parar la tabla y trandofrmarla en data.frame
cts.norm=t(cts.trans.norm) %>% 
  as.data.frame()



## -------------------------------------------------------------------------
topgenes.n=5 # cantidad de genes a considerar


## -------------------------------------------------------------------------

genes=rownames(res2)[1:topgenes.n] # nombres de estos genes

genes.exp=cts.norm[genes,] %>% # extracción de estos genes de la tabla de conteos
  t() # verticalización para poder juntar con la tabla del diseño experimental

exp.table=cbind(samples, genes.exp) # unión con la tabla de diseño

exp.table.long=pivot_longer(exp.table, cols = -c(1:5), 
             names_to = 'gene', values_to = 'exp') # verticalización




## -------------------------------------------------------------------------
data.plot= as.data.frame(exp.table.long) %>% # cambiamo a CLASS data.frame para trabajar con filter()
         filter(Treatment %in% trat) %>% # filtramos por tratamiento (sólo ctr)
         filter(Tissue %in% tejidos) # filtramos por tejido (sólo CHR)) 



## ---- fig.height=5, fig.width=3, dev='svg', dpi=300, fig.align='center'----

ggplot(data= data.plot) +  # datos en los que ggplot buscará lo que le ordenemos en la progrmación del gráfico
  
  # generamos un elemento (geom) de puntos
  geom_jitter(aes(x=Sex, y=exp, # la posición de cada punto codificará el nivel de expresión en función del sexo
                  fill=Sex), # el relleno de los puntos también codificará el sexo
              width = 0.05,  # ancho de la dispersión horizontal aleatoria (para que no se solapen tanto)
              height = 0,
             pch=21, # estilo de punto
             size=3, # tamaño de los puntos
             alpha=0.5)+ # transparencia (0 - 1)
  
  scale_fill_manual(values = c(M='slateblue2', F='hotpink1'))+ # colores de cada sexo
  
  facet_grid(gene~., # desdoblamos el gráfico verticalmente, uno para cada gen
             labeller = labeller(gene = paste0(
               # título de cada panel
             genes, '
LFC = ', res$log2FoldChange[rownames(res) %in% genes] %>% round(2)) %>% 
               `names<-`(genes))
               ,
             scales = 'free')+ # cada panel puede tener su propia escala
  
  labs(x=NULL, y='Expresión')+ # títulos de ejes
  
  theme_bw() # estética general del gráfico




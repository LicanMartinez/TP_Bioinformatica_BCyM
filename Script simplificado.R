#### ANALIZANDO DATOS DE RNASEQ ----------------------------------
#
#    SCRIPT PARA GUIA DE Tp 4 - Biolog?a Celular y Molecular
#    Centro Regional Universitario Bariloche
#    Universidad Nacional del Comahue
#
#    Generado por Eduardo E. Zattara
#    a partir de la vi?eta escrita por
#    Michael I. Love, Simon Anders, Vladislav Kim and Wolfgang Huber. 2019.
#    RNA-seq workflow: gene-level exploratory analysis and differential expression. 
#    http://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
#

#### Script Version 1.1 -----

#### 3.3 - Instalando los paquetes requeridos -----------------------
# Paquetes a instalar
#install.packages(c("dplyr","tidyr","ggplot2","pheatmap",
#"RColorBrewer","PoiClaClu","glmpca","ggbeeswarm"))



# glmpca requiere R 3.6+

### SI LA VERSION DE R < 3.5
# Instalar paquetes de Bioconductor

# source("https://bioconductor.org/biocManager.R")
# BiocInstaller::biocManager(c("DESeq2", "vsn","apeglm","genefilter","IHW","edgeR"))


### SI LA VERSION DE R >= 3.5

install.packages("BiocManager")

#BiocManager::install(version = "3.12", force = T)

#BiocManager::install(c("DESeq2", "vsn","apeglm","genefilter","IHW", 'edgeR'), force = T)

#Invocar bibliotecas 
l=c("dplyr", 
    "tidyr",
    "ggplot2",
    "pheatmap",
    "RColorBrewer",
    "PoiClaClu",
    "glmpca",
    "ggbeeswarm", 
    "tidyr", 
    "DESeq2")

lapply(l, require, character.only = TRUE)


####3.4 - Crear un objeto DESeqDataSet a partir de la matriz de conteos y la tabla de muestras---- 

#!!! Es fundamental que las columas de la matriz de conteos y las filas de las
# muestras est?n en el mismo orden.

#Leer tabla de conteos
dsxCts <- read.delim("GSE87788_OnthophagusCounts.tsv", row.names = "gene")

#Leer la tabla de informaci?n sobre las muestras, y seleccionar variables de inter?s
dsxSamples <- read.delim("dsxSamples.tsv", sep="\t", row.names=1)
dsxColdata <- dsxSamples[,c("Tissue","Treatment","Sex","individual")]

#Verificar que todas las muestras de la tabla est?n presentes con el mismo nombre en la matriz
all(rownames(dsxColdata) %in% colnames(dsxCts)) # Debe devolver "TRUE"
all(rownames(dsxColdata) == colnames(dsxCts)) # Debe devolver "TRUE"

#Creamos el objeto DESeqDataSet
dds_dsx_all <- DESeqDataSetFromMatrix(countData = dsxCts,
                                       colData = dsxColdata,
                                       design = ~ Tissue + Sex + Treatment)


#### 4.1- Prefiltrando el conjunto de datos ---------------
#Filtramos para remover genes con conteos muy bajos (<1)
nrow(dds_dsx_all)
keep <- rowSums(counts(dds_dsx_all)) > 1

#Filtramos para remover genes con conteos muy bajos (<10)
keep <- rowSums(counts(dds_dsx_all)) >= 10
dds_dsx_work <- dds_dsx_all[keep,]
nrow(dds_dsx_work)

#### 4.2 - La transformaci?n estabilizadora de la varianza -----

# Usando vst para estabilizar la varianza 
# el nuevo objeto vsd tiene informaci?n equivalente a dds_dsx_work pero con los valores 
# de conteos procesados para que tengan propiedades estad?sticas m?s manejables
vsd <- vst(dds_dsx_work, blind = FALSE)

#### 4.3 - Distancias entre muestras --------

## Euclideas (con datos vsd)
# Transponemos la tabla de conteos transformados porque dist espera muestras en filas
sampleDists <- dist(t(assay(vsd)))

library("pheatmap")
library("RColorBrewer")

#Preparamos la matriz de datos, le asignamos nombres a cada muestra, y corremos el heatmap
# la matriz tiene en cada celda el valor de distancia entre el par de muestras fila-clumna

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste((vsd)$Tissue, (vsd)$Sex,  (vsd)$Treatment, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, fontsize = 7)

## De Poisson (con conteos crudos)
library("PoiClaClu")
#Primero transponemos los conteos dds_dsx_work
poisd <- PoissonDistance(t(counts(dds_dsx_work)))

#Trazamos el mapa de calor en una figura a continuaci?n
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste((vsd)$Tissue, (vsd)$Sex,  (vsd)$Treatment, sep = " - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors, fontsize = 7)



#### 4.4 - Visualizando distancias entre muestras usando Componentes Principales ------
#Visualizacion usando plotPCA 
plotPCA(vsd, intgroup = c("Tissue", "Sex", "Treatment"))

# (Opcional) Visualizacion usando ggplot (m?s complicado pero m?s personalizable)
#Volcamos los datos de plotPCA a un objeto
pcaData <- plotPCA(vsd, intgroup = c("Tissue", "Sex", "Treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
SexTreatment <-with(pcaData, paste(Sex,Treatment,sep="-"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = SexTreatment, shape = Tissue)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")

#### 5 - Analisis de Expresion Diferencial ------------

####5.1 - Correr analisis general (compara tratamiento control vs dsx RNAi)----
dds_dsx_work <- estimateSizeFactors(dds_dsx_work)
dds_dsx_work <- DESeq(dds_dsx_work)

####5.2 - Construyendo la tabla de resultados-----

res <- results(dds_dsx_work, contrast=c("Treatment","dsx","ctr"))

#Resumen del testeo estadistico
summary(res)

#Ajustando los umbrales de p-valores
res.05 <- results(dds_dsx_work, alpha = 0.05)
table(res.05$padj < 0.05)

# y de log-fold-change
resLFC1 <- results(dds_dsx_work, lfcThreshold=1)
table(resLFC1$padj < 0.1)

#### 5.3 - Otras comparaciones -----

# Contrastes espec?ficos

summary(results(dds_dsx_work, contrast = c("Tissue", "BRN", "CHE")))

#### 5.4 - Testeos m?ltiples ----

#Vamos a seleccionar solo aquellos genes con un p ajustado <0.1
resSig <- subset(res, padj < 0.1)

# Ordemamos por log-fold-chanche de menor a mayor
# (con head vemos solo las primears filas)
filas=5 # elegimos la cantidad de filas

head(resSig[ order(resSig$log2FoldChange), ], filas)

# y de mayor a menor
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ], filas)

#### 6 - Graficando los resultados -------------

#### 6.1 - Gr?fico de conteos -------
#Calculamos resultados para el contraste entre cerebro y epitelio cef?lico
brn_che_res <- results(dds_dsx_work, contrast = c("Tissue", "BRN", "CHE"))
#Buscamos el nombre del gen con valor m'as bajo de p ajustado 
topGene <- rownames(brn_che_res)[which.min(brn_che_res$padj)]
#Usamos plotCounts para mostrar los conteos normalizados para ese gen
plotCounts(dds_dsx_work, gene = topGene, intgroup=c("Tissue"), pch=16)

#Haciendo graficos personalizados usando ggplot2
library("ggbeeswarm")
geneCounts <- plotCounts(dds_dsx_work, gene = topGene, intgroup = c("Tissue","Sex","Treatment"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = Tissue, y = count, color = Sex, shape = Treatment)) +
  geom_beeswarm (cex = 2.5, alpha=0.7, size=3)+scale_y_log10()+
  ggtitle(topGene)

#### 6.2 - Gr?ficos MA -----

#Primero, reducimos los log2fold changes

#library("apeglm")
#res <- lfcShrink(dds_dsx_work, coef="Treatment_dsx_vs_ctr", type="apeglm")
#plotMA(res, ylim = c(-2, 2))

#Equivalente sin usar lfcShrink
res.noshr <- results(dds_dsx_work, name="Treatment_dsx_vs_ctr")
plotMA(res.noshr, ylim = c(-2, 2))

#Etiquetando genes particulares en el gr?fico MA
plotMA(res.noshr, ylim = c(-2,2))
topGene <- rownames(res.noshr)[which.min(res.noshr$padj)]
with(res.noshr[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

#Histograma de p values
hist(res$pvalue[res$baseMean > 1], breaks = 0:50/50,
     col = "grey50", border = "white")

####6.3 - Clustering de genes -----
library("genefilter")

#Seleccionamos los 20 genes m?s variables
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)

#y armamos un heatmap
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Tissue","Sex","Treatment")])
pheatmap(mat, annotation_col = anno, fontsize_col = 6)

#### 6.4 - Filtrado independiente ------

# Usando la tabla resLFC1 generada previamente
qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(resLFC1$pvalue, bins, function(p) mean(p < .05, na.rm = TRUE))
barplot(fractionSig, xlab = "mean normalized count", ylab = "fraction of small p values")

#### 6.5 - Ponderaci?n de hip?tesis independiente
library("IHW")
res.ihw <- results(dds_dsx_work, filterFun=ihw)

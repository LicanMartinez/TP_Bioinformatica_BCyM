#### ANALIZANDO DATOS DE RNASEQ ----------------------------------
#
#    SCRIPT PARA GUIA DE Tp 4 - Biología Celular y Molecular
#    Centro Regional Universitario Bariloche
#    Universidad Nacional del Comahue
#
#    Generado por Eduardo E. Zattara
#    a partir de la viñeta escrita por
#    Michael I. Love, Simon Anders, Vladislav Kim and Wolfgang Huber. 2019.
#    RNA-seq workflow: gene-level exploratory analysis and differential expression. 
#    http://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
#

#### Script Version 1.1 -----

#### 3.3 - Instalando los paquetes requeridos -----------------------
# Paquetes a instalar
install.packages(c("dplyr", 
                   "tidyr",
                 "ggplot2",
                 "pheatmap",
                 "RColorBrewer",
                 "PoiClaClu",
                 "glmpca",
                 "ggbeeswarm"))

# glmpca requiere R 3.6+

### SI LA VERSION DE R < 3.5
# Instalar paquetes de Bioconductor
#source("https://bioconductor.org/biocManager.R")
#BiocInstaller::biocManager(c("DESeq2", "vsn","apeglm",
#                             "genefilter","IHW","edgeR"))

### SI LA VERSION DE R >= 3.5
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")

BiocManager::install(c("DESeq2", "vsn","apeglm","genefilter","IHW"))
BiocManager::install(c("DESeq2", "vsn","apeglm","genefilter","IHW"))

#Invocar bibliotecas 
library("GenomicFeatures")
library(tidyr)
library("DESeq2")

####3.4 - Crear un objeto DESeqDataSet a partir de la matriz de conteos y la tabla de muestras---- 

#!!! Es fundamental que las columas de la matriz de conteos y las filas de las
# muestras estén en el mismo orden.

#Leer tabla de conteos
dsxCts <- read.delim("GSE87788_OnthophagusCounts.tsv", row.names = "gene")

#Leer la tabla de información sobre las muestras, y seleccionar variables de interés
dsxSamples <- read.delim("dsxSamples.tsv", sep="\t", row.names=1)
dsxColdata <- dsxSamples[,c("Tissue","Treatment","Sex","individual")]

#Verificar que todas las muestras de la tabla estén presentes con el mismo nombre en la matriz
all(rownames(dsxColdata) %in% colnames(dsxCts)) # Debe devolver "TRUE"
all(rownames(dsxColdata) == colnames(dsxCts)) # Debe devolver "TRUE"

#Creamos el objeto DESeqDataSet
dds_dsx_all <- DESeqDataSetFromMatrix(countData = dsxCts,
                                       colData = dsxColdata,
                                       design = ~ Tissue + Sex + Treatment)


#### 4.1- Prefiltrando el conjunto de datos ---------------
#Filtramos para remover genes con conteos muy bajos (<1)
nrow(dds_dsx_all)
keep <- rowSums(counts(dds_dsx_all)) > 10

#Filtramos para remover genes con conteos muy bajos (<10)
keep <- rowSums(counts(dds_dsx_all)) >= 10
dds_dsx_work <- dds_dsx_all[keep,]
nrow(dds_dsx_work)

#### 4.2 - La transformación estabilizadora de la varianza y el rlog -----

#Viendo por qué es importante estabilizar la varianza 
lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library("vsn")
meanSdPlot(cts, ranks = FALSE)

#Y para los recuentos transformados por logaritmo:
log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)

# Usando vst para estabilizar la varianza 
vsd <- vst(dds_dsx_work, blind = FALSE)
head(assay(vsd)[,1:4], 3) 
colData(vsd)

#### Usando rlog para estabilizar la varianza 
#rlog() may take a long time with 50 or more samples,
#vst() is a much faster transformation!
rld <- rlog(dds_dsx_work, blind = FALSE)
head(assay(rld)[,1:4], 3)

#### Y ahora comparamos ambos enfoques 
library("dplyr")
library("ggplot2")

dds_dsx_work <- estimateSizeFactors(dds_dsx_work)
df <- bind_rows(
  as_data_frame(log2(counts(dds_dsx_work, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation) 

#### 4.3 - Distancias entre muestras --------

#### 4.3.1 - Distancia euclideana entre muestras----------------
# Transponemos la tabla de conteos transformados porque dist espera muestras en filas
sampleDists <- dist(t(assay(vsd)))

library("pheatmap")
library("RColorBrewer")

#Preparamos la matriz de datos, le asignamos nombres a cada muestra, y corremos el heatmap
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste((vsd)$Tissue, (vsd)$Sex,  (vsd)$Treatment, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, fontsize = 7)

#### 4.3.2 - Calculando distancias usando distancias de Poisson
library("PoiClaClu")
#Primero transponemos los conteos dds_dsx_work
poisd <- PoissonDistance(t(counts(dds_dsx_work)))

#Trazamos el mapa de calor en una figura a continuación
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste((vsd)$Tissue, (vsd)$Sex,  (vsd)$Treatment, sep = " - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors, fontsize = 7)

#### 4.4 - Visualizando distancias entre muestras usando Componentes Principales ------
#Visualizacion usando plotPCA 
plotPCA(vsd, intgroup = c("Tissue", 'Treatment'))

#Graficacion usando ggplot
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

#### 4.5 - Visualizando distancias entre muestras usando Componentes Principales generalizados------
library("glmpca") # requiere R versión 3.6+
gpca <- glmpca(counts(dds_dsx_work), L=2)
gpca.dat <- gpca$factors
gpca.dat$Sex <- dds_dsx_work$Sex
gpca.dat$Tissue <- dds_dsx_work$Tissue
gpca.dat$Treament <- dds_dsx_work$Treatment
gpca.dat$SexTreatment <- paste(dds_dsx_work$Sex,dds_dsx_work$Treatment,sep="-")
ggplot(gpca.dat, aes(x = dim1, y = dim2, color = SexTreatment, shape = Tissue)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

#### 4.6 - Visualizando distancias entre muestras usando MultiDimensional Scaling------
# Usando distancias euclidianas
mds <- as.data.frame(colData(vsd))  %>% cbind(cmdscale(sampleDistMatrix))
SexTreatment <-with(mds, paste(Sex,Treatment,sep="-"))
ggplot(mds, aes(x = `1`, y = `2`, color = SexTreatment, shape = Tissue)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

# Usando distancias de Poisson
mdsPois <- as.data.frame(colData(dds_dsx_work)) %>% cbind(cmdscale(samplePoisDistMatrix))
SexTreatment <-with(mdsPois, paste(Sex,Treatment,sep="-"))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = SexTreatment, shape = Tissue)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")

#### 5 - Analisis de Expresion Diferencial ------------

####5.1 - Correr analisis general (compara tratamiento control vs dsx RNAi)----
dds_dsx_work <- DESeq(dds_dsx_work)

####5.2 - Construyendo la tabla de resultados-----

#Llamar a la funcion de resultados sin parametros adicionales
res <- results(dds_dsx_work)

#Ahora especificando un contraste
res <- results(dds_dsx_work, contrast=c("Treatment","dsx","ctr"))

#Viendo los metadatos de los resultados
mcols(res, use.names = TRUE)

#Resumen del testeo estadistico
summary(res)

#Ajustando los umbrales
res.05 <- results(dds_dsx_work, alpha = 0.05)
table(res.05$padj < 0.05)

resLFC1 <- results(dds_dsx_work, lfcThreshold=1)
table(resLFC1$padj < 0.1)

#### 5.3 - Otras comparaciones -----

# Contrastes específicos
results(dds_dsx_work, contrast = c("Tissue", "BRN", "CHE"))
summary(results(dds_dsx_work, contrast = c("Tissue", "BRN", "CHE")))

#### 5.4 - Testeos múltiples ----
#Explorando la importancia del FDR (false discovery ratio)
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))
sum(res$padj < 0.1, na.rm=TRUE)

#Vamos a seleccionar solo aquellos genes con un p ajustado <0.1 y ordenarlos
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

#### 5.5 - Generando una tabla de conteos normalizados -----
library(edgeR)
#Primero descartamos todas las filas con conteo total de 0
dsxCts_nozeros <- dsxCts[which(rowSums(dsxCts)>0),]

#Generamos la tabla con normalizacion TMM (normaliza por  tamaño de biblioteca, pero no por largo de gen)
rnaseqMatrix = as.matrix(dsxCts_nozeros)
exp_study = DGEList(counts=rnaseqMatrix, group=factor(colnames(rnaseqMatrix)))
exp_study = calcNormFactors(exp_study)
exp_study$samples$eff.lib.size = exp_study$samples$lib.size * exp_study$samples$norm.factors
dsx_TMM_matrix <- cpm(rnaseqMatrix)

#Exportamos la tabla como archivo de texto delimitado por tabulador
write.table(dsx_TMM_matrix, file="dsx.TMM_table.tsv",quote = F,sep = "\t",row.names = T, col.names = T)


#### 6 - Graficando los resultados -------------

#### 6.1 - Gráfico de conteos -------
#Calculamos resultados para el contraste entre cerebro y epitelio cefálico
brn_che_res <- results(dds_dsx_work, contrast = c("Tissue", "BRN", "CHE"))
#Buscamos el nombre del gen con valor m'as bajo de p ajustado 
topGene <- rownames(brn_che_res)[which.min(brn_che_res$padj)]
#Usamos plotCounts para mostrar los conteos normalizados para ese gen
plotCounts(dds_dsx_work, gene = topGene, intgroup=c("Tissue"), pch=16)

#Haciendo graficos personalizados usando ggplot2
library("ggbeeswarm")
geneCounts <- plotCounts(dds_dsx_work, gene = topGene, intgroup = c("Tissue","Sex","Treatment"),returnData = TRUE)
ggplot(geneCounts, aes(x = Tissue, y = count, color = Sex, shape = Treatment)) +
  geom_beeswarm (cex = 3)+   scale_y_log10() 

#### 6.2 - Gráficos MA -----

#Primero, reducimos los log2fold changes

library("apeglm")
resultsNames(dds_dsx_work)
## [1] "Intercept"               "cell_N061011_vs_N052611"
## [3] "cell_N080611_vs_N052611" "cell_N61311_vs_N052611" 
## [5] "dex_trt_vs_untrt"
res <- lfcShrink(dds_dsx_work, coef="Treatment_dsx_vs_ctr", type="apeglm")
plotMA(res, ylim = c(-2, 2))

#Equivalente sin usar lfcShrink
res.noshr <- results(dds_dsx_work, name="Treatment_dsx_vs_ctr")
plotMA(res.noshr, ylim = c(-2, 2))

#Etiquetando genes particulares en el gráfico MA
plotMA(res, ylim = c(-2,2))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

#Histograma de p values
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

####6.3 - Clustering de genes -----
library("genefilter")

#Seleccionamos los 20 genes más variables
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

#### 6.5 - Ponderación de hipótesis independiente
library("IHW")
res.ihw <- results(dds_dsx_work, filterFun=ihw)

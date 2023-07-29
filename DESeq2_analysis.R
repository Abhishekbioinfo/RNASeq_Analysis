#!/usr/bin/env Rscript

rm(list=ls()) 
suppressMessages(library(vidger))
suppressMessages(library(DESeq2))
#suppressMessages(library(xlsx))
suppressMessages(library(DT))
suppressMessages(library(ggplot2))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(regionReport))
suppressMessages(library(DEFormats))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
suppressMessages(library(dplyr))
suppressMessages(library(colorspace))
suppressMessages(library(optparse))
suppressMessages(library(scales))
suppressMessages(library(readr))
suppressMessages(library(DelayedArray))
#suppressMessages(library(apeglm))
suppressMessages(library(gplots))

library(reshape)
library(ggrepel)
library(DEGreport)
#library(RVA)

#BiocManager::install("GSVAdata")
#BiocManager::install("ComplexHeatmap")
#BiocManager::install("rWikiPathways")


                                 # to run the script in command lines

# options list with associated default value.
option_list <- list( 
make_option(c("-n", "--projectName"),
			default=basename(getwd()),
			dest="projectName",
			help="name of the project used for storing images and tables [default: name of the current directory]."),

make_option(c("-d", "--design"),
			default="design.txt",
			dest="design",
			help="path to the design/target file [default: %default]."),

make_option(c("-q", "--featureQuant"),
			dest="featureQuant",
			help="path to the featureCount Quantification file [default: %default]."),

make_option(c("-g", "--geneFeatureTable"),
      dest="geneFeatureTable",
      help="path to the Feature Table file containg the following headers  [default: %default]."),

	
make_option(c("-v", "--varInt"),
			default="group",
			dest="varInt", 
			help="factor of interest [default: %default]"),

make_option(c("-c", "--condControl"),
			default="control",
			dest="condControl",
			help="control biological condition [default: %default]"),
			
make_option(c("-t", "--condTreated"),
			default="treated",
			dest="condTreated",
			help="treated biological condition [default: %default]"),

make_option(c("-b", "--batch"),
			default=NULL,
			dest="batch",
			help="blocking factor [default: %default] or \"batch\" for example"),

make_option(c("-f", "--fitType"),
			default="parametric",
			dest="fitType", 
			help="mean-variance relationship: [default: %default],local or mean"),

make_option(c("-a", "--alpha"),
			default=0.05,
			dest="alpha", 
			help="threshold of statistical significance [default: %default]"),

make_option(c("-p", "--pAdjustMethod"),
			default="BH",
			dest="pAdjustMethod", 
			help="p-value adjustment method: \"BH\" or \"BY\" [default: %default]"),


make_option(c("-l", "--locfunc"),
			default="median",
			dest="locfunc", 
			help="median or shorth to estimate the size factors [default: %default]")

)

# now parse the command line to check which option is given and get associated values
parser <- OptionParser(usage="usage: %prog [options]",
					   option_list=option_list, 
					   description="Compare two or more biological conditions in a RNA-Seq framework with DESeq2.",
					   epilogue="Computational Genome Biology Lab, CSIR-IICB, Kolkata 700032")
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options



#Check mandetory inputs 
if ( is.null(opt$design) ) {
  stop("--sample groupfile file / target file must be provided. See script usage (--help)")
}

if ( is.null(opt$featureQuant) ) {
  stop("--path to the featureCounts Countfile file must be provided. See script usage (--help)")
}


if ( is.null(opt$condControl) ) {
  stop("--reference / control biological condition name must be provided. See script usage (--help)")
}


if ( is.null(opt$condTreated) ) {
  stop("--treatment biological condition name must be provided. See script usage (--help)")
}



# get options and arguments
workDir <- getwd()
projectName <- opt$projectName 
reportName <-opt$reportName                          # name of the report
design <- opt$design   		       # path to the design/target file   
featureQuant <- opt$featureQuant                     # path to the featureCounts file
varInt <- opt$varInt                                 # factor of interest
condControl <- opt$condControl   		       # reference biological condition
condTreated <- opt$condTreated                       # treated biological condition
batch <- opt$batch                                   # blocking factor: NULL (default) or "batch" for example
fitType <- opt$fitType                               # mean-variance relationship: "parametric" (default), "local" or "mean"
alpha <- as.numeric(opt$alpha)                       # threshold of statistical significance
pAdjustMethod <- opt$pAdjustMethod                   # p-value adjustment method: "BH" (default) or "BY"
locfunc <- opt$locfunc                               # "median" (default) or "shorth" to estimate the size factors
				

print(paste("workDir", workDir))
print(paste("projectName", projectName))
print(paste("design", design))
print(paste("featureQuant", featureQuant))
print(paste("varInt", varInt))
print(paste("condControl", condControl))
print(paste("condTreated", condTreated))
print(paste("batch", batch))
print(paste("fitType", fitType))
print(paste("alpha", alpha))
print(paste("pAdjustMethod", pAdjustMethod))
print(paste("locfunc", locfunc))

################################################################################
###                             running script                               ###
################################################################################

#dir.create(projectName, showWarnings = FALSE, recursive = TRUE)
#dir.create(projectName, "tables", showWarnings = FALSE, recursive = TRUE)
#dir.create(projectName, "figures", showWarnings = FALSE, recursive = TRUE)

dir.create(file.path(projectName,"tables"), recursive = TRUE)
dir.create(file.path(projectName,"figures"), recursive = TRUE)


#Check design File
checkdesignFile <- function(design, varInt, condControl, condTreated, batch){
  target <- read.table(design, header=TRUE, sep="\t", na.strings="")
  if (!I(varInt %in% names(target))) stop(paste("The factor of interest", varInt, "is not in the target file"))
  if (!is.null(batch) && !I(batch %in% names(target))) stop(paste("The batch effect", batch, "is not in the target file")) 
  target[,varInt] <- as.factor(target[,varInt])
  if (!I(condControl %in% as.character(target[,varInt]))) stop(paste("The reference level", condControl, "is not a level of the factor of interest"))
  if (!I(condTreated %in% as.character(target[,varInt]))) stop(paste("The reference level", condTreated, "is not a level of the factor of interest"))
  target[,varInt] <- relevel(target[,varInt],ref=condControl)
  target <- target[order(target[,varInt]),]
  rownames(target) <- as.character(target[,1])
  # check if varInt contains replicates
  if (min(table(target[,varInt]))<2) stop(paste("The factor of interest", varInt, "has a level without replicates"))
  # check if NA in the target
  if (any(is.na(cbind(target[,c(varInt, batch)], target[,1:2])))) stop("NA are present in the target file")
  # warning message if batch is numeric
  if (!is.null(batch) && is.numeric(target[,batch])) warning(paste("The", batch, "variable is numeric. Use factor() or rename the levels with letters to convert it into a factor"))
  if (any(grepl("[[:punct:]]", as.character(target[,varInt])))) stop(paste("The", varInt, "variable contains punctuation characters, please remove them"))
  cat("Target file:\n")
  print(target)
  return(target)
}


target <- checkdesignFile(design=design, varInt=varInt, condControl=condControl, condTreated=condTreated, batch=batch)

				   
# loading target file
target <- checkdesignFile(design=design, varInt=varInt, condControl=condControl, condTreated=condTreated, batch=batch)
rownames(target) <- as.character(target[,1])

countData <- read.table(featureQuant, header=TRUE, sep="\t", row.names=1)
transriptData <- countData[,1:5]
countData <- countData[ ,6:ncol(countData)]

colnames(countData) <- gsub("\\.[sb]am$", "", colnames(countData))
colnames(countData) <- lapply(colnames(countData), function(x) sapply(strsplit(x, "\\."), tail, 1))
colnames(countData)
countData <- as.matrix(countData)

all(rownames(target) %in% colnames(countData))

countData <- countData[, rownames(target)]
all(rownames(target) == colnames(countData))

write.table(countData, file=paste0(projectName,"/tables/","counts.tsv"), sep="\t", dec=".", quote=FALSE,col.names=NA)



run.DESeq2 <- function(counts, target, varInt, batch=NULL,
                       locfunc="median", fitType="parametric", pAdjustMethod="BH",alpha=0.05, ...){
  # building dds object
  dds <- DESeqDataSetFromMatrix(countData=countData, colData=target, 
                                design=formula(paste("~", ifelse(!is.null(batch), paste(batch,"+"), ""), varInt)))
  cat("Design of the statistical model:\n")
  cat(paste(as.character(design(dds)),collapse=" "),"\n")           
  
  # normalization
  dds <- estimateSizeFactors(dds,locfunc=eval(as.name(locfunc)))
  cat("\nNormalization factors:\n")
  print(sizeFactors(dds))
  
  # estimating dispersions
  dds <- estimateDispersions(dds, fitType=fitType)
  
  # statistical testing: perform all the comparisons between the levels of varInt
  dds <- nbinomWaldTest(dds, ...)

  results <- list()
  for (comp in combn(nlevels(colData(dds)[,varInt]), 2, simplify=FALSE)){
    levelRef <- levels(colData(dds)[,varInt])[comp[1]]
    levelTest <- levels(colData(dds)[,varInt])[comp[2]]
    results[[paste0(levelTest,"_vs_",levelRef)]] <- results(dds, contrast=c(varInt, levelTest, levelRef),
                                                            pAdjustMethod=pAdjustMethod, alpha=alpha)
    cat(paste("Comparison", levelTest, "vs", levelRef, "done\n"))
  }
  
  return(list(dds=dds,results=results,sf=sizeFactors(dds)))
}

out.DESeq2 <- run.DESeq2(counts=countData, target=target, varInt=varInt, batch=batch,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         alpha=alpha)

dds <- out.DESeq2$dds

mcols(dds)$basepairs <- transriptData[,5]
fpkm <- data.frame(GeneID=row.names(countData),length=transriptData[,5],fpkm(dds,robust=T)) 

write.table(fpkm, file=paste0(projectName,"/tables/",condTreated,"_vs_",condControl,"_FPKM.tsv"), sep="\t", dec=".", quote=FALSE, row.names = FALSE)



###########################
exportResults.DESeq2 <- function(out.DESeq2, group, alpha=0.05, export=TRUE){
  
  dds <- out.DESeq2$dds
  results <- out.DESeq2$results
  
  # comptages bruts et normalis?s
  counts <- data.frame(Id=rownames(counts(dds)), counts(dds), round(counts(dds, normalized=TRUE)))
  colnames(counts) <- c("ID", colnames(counts(dds)), paste0("norm.", colnames(counts(dds))))
  # baseMean avec identifiant
  bm <- data.frame(ID=rownames(results[[1]]),baseMean=round(results[[1]][,"baseMean"],2))
  # merge des info, comptages et baseMean selon l'Id
  base <- merge(counts, bm, by="ID", all=TRUE)
  tmp <- base[,paste("norm", colnames(counts(dds)), sep=".")]
  for (cond in levels(group)){
    base[,cond] <- round(apply(as.data.frame(tmp[,group==cond]),1,mean),0)
  }
  
  complete <- list()
  for (name in names(results)){
    complete.name <- base

    # ajout d'elements depuis results
    res.name <- data.frame(ID=rownames(results[[name]]),
                           FoldChange=round(2^(results[[name]][,"log2FoldChange"]), 3),
                           log2FoldChange=round(results[[name]][,"log2FoldChange"], 3),
                           stat=round(results[[name]][,"stat"], 3),
                           pvalue=results[[name]][,"pvalue"],
                           padj=results[[name]][,"padj"])
    complete.name <- merge(complete.name, res.name, by="ID", all=TRUE)
    # ajout d'elements depuis mcols(dds)
    mcols.add <- data.frame(ID=rownames(counts(dds)),dispGeneEst=round(mcols(dds)$dispGeneEst,4),
                            dispFit=round(mcols(dds)$dispFit,4),dispMAP=round(mcols(dds)$dispMAP,4),
                            dispersion=round(mcols(dds)$dispersion,4),betaConv=mcols(dds)$betaConv,
                            maxCooks=round(mcols(dds)$maxCooks,4))
    complete.name <- merge(complete.name, mcols.add, by="ID", all=TRUE)
    complete[[name]] <- complete.name
    
    if (export){
      # s?lection des up et down
      # write.table(complete.name, file=paste0("tables/",name,".complete.txt"), sep="\t", row.names=FALSE, dec=".", quote=FALSE)
      write.table(complete.name, file=paste0(projectName,"/tables/",condTreated,"_vs_",condControl,"_complete_DESeq2_result.tsv"), sep="\t", dec=".", quote=FALSE, row.names = FALSE)
     
    }
  }

  return(complete)
}


###########################
exportResults.DESeq2(out.DESeq2, group=unique(target$group), alpha=alpha)


#res <- results(dds, name=paste0(condTreated,"_vs_",condControl))


coldata <- colData(dds)
print("#### COLDATA ######")
coldata
intgroup<- colnames(coldata[2])


print("####  ######")
intgroup

dds.rld.trans <- rlog(dds, blind=FALSE)


rld <- assay(rlog(dds, blind=FALSE))

vsd <- vst(dds, blind = FALSE)



#exportResults.DESeq2(out.DESeq2, group=unique(target$group), alpha=alpha)

#save.image(file=paste0(projectName,"/", projectName,".RData"))
#write.table(rld, file=paste0(projectName,"/",condTreated,"_vs_",condControl,"_RLD.tsv"), sep="\t", dec=".", quote=FALSE,col.names=NA)
#write.table(fpm, file=paste0(projectName,"/",condTreated,"_vs_",condControl,"_FPM.tsv"), sep="\t", dec=".", quote=FALSE,col.names=NA)


#Plots

## For ggplot

res<- results(dds, contrast=c(varInt, condTreated, condControl), pAdjustMethod=pAdjustMethod, alpha=alpha)

res_sig <- subset(res, padj<=.05)
res_lfc_1 <- subset(res_sig, abs(log2FoldChange) > 1) 
res_lfc_2 <- subset(res_sig, abs(log2FoldChange) > 2)
lfc_2 <- subset(res, abs(log2FoldChange) > 2)

res.df <- as.data.frame(res)
ord <- order(res.df$padj, decreasing = FALSE)
res.df <- res.df[ord, ]
features <- rownames(res.df)
res.df <- cbind(data.frame(Feature = features), res.df)
rownames(res.df) <- NULL


normalized_counts <- counts(dds, normalized=T)
norm_OEsig <- normalized_counts[rownames(res_lfc_1),]

annotation <- data.frame(condition=target[, varInt], 
                     row.names=rownames(target))

### Set a color palette
heat.colors <- brewer.pal(6, "Blues") 
#YlOrRd
### Run pheatmap
gene_heatmap<-pheatmap(norm_OEsig, color = heat.colors, cluster_rows = T, show_rownames=F,
annotation= annotation, border_color=NA, fontsize = 10, scale="row",
     fontsize_row = 10, height=20)

ggsave(file=paste0(projectName,"/figures/","gene_heatmap.png"), gene_heatmap, width = 18, dpi = 300, units = "cm", device='png')



#Sample Distance
sampleDists <- dist(t(assay(vsd)))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)
sampleDistMatrix
rownames(sampleDistMatrix) <- paste(  vsd$group, vsd$samples,sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#png("heatmap-samples.png", w=1000, h=1000, pointsize=20)

sample_dist <- pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors, fontsize = 12)
ggsave(file=paste0(projectName,"/figures/","SAMPLE_Distance.png"), sample_dist, width = 18, dpi = 300, units = "cm", device='png')
#dev.off()

#ggsave(file="Sample_to_Sample_Distance.png", sample_distance_plot, width = 15, dpi = 300, units = "cm", device='png')


## Transform count data
rld <- tryCatch(rlog(dds), error = function(e) { rlog(dds, fitType = 'mean') })
## Perform PCA analysis and make plot
#png("PCA.png", w=1000, h=1000, pointsize=20)


nudge <- position_nudge(y = 2)
PCA <-plotPCA(rld, intgroup = intgroup) + geom_text(aes(label = name), position = nudge) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggsave(file=paste0(projectName,"/figures/","PCA.png"), PCA, width = 18, dpi = 300, units = "cm", device='png')


#Gene Cluster

topVarGenes <- head(order(-rowVars(assay(rld))),50)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("samples",varInt)])
print("")
print("Generating Gene Heatmap")
print("")
geneheatmap<-pheatmap(mat,annotation_col=df, fontsize = 6)
ggsave(file=paste0(projectName,"/figures/","gene_heatmap_plot.png"), geneheatmap, width = 18, dpi = 300, units = "cm", device='png')

#Volcanoplot

print("")
print("Generating Enhance Volcano Plot")
print("")
#png("VOLCANO.png", w=1000, h=1000, pointsize=15)
evp <- EnhancedVolcano(res,
lab =  rownames(res),
x = 'log2FoldChange',
y = 'padj',
xlab = bquote(~Log[2]~ 'fold change'),
pCutoff = 0.05,
FCcutoff = 2.0,
pointSize = 1.0,
labSize = 3.0,
colAlpha = 1,
title = NULL,
subtitle = NULL,
legendPosition = 'top',
legendLabSize = 12,
legendIconSize = 4.0,
drawConnectors = TRUE,
widthConnectors = 0.75) + xlim(-15,15)
ggsave(file=paste0(projectName,"/figures/","Enhanced_volcano_plot.png"), evp, width = 18, dpi = 300, units = "cm", device='png')


#MAPlot
png(file=paste0(projectName,"/figures/","DESeq2-MA-Plot.png"), w=1000, h=1000, pointsize=20)
plotMA(res, ylim = c(-5,5))
#ggsave(file=paste0(projectName,"/figures/","DESeq-MA.png"), deseqMA, width = 18, dpi = 300, units = "cm", device='png')
dev.off()


vmaplot <- vsMAPlot(
    x = condControl, y = condTreated, 
    data = dds, d.factor = varInt, type = 'deseq', 
    padj = 0.05, y.lim = NULL, lfc = 1, title = TRUE, 
    legend = TRUE, grid = TRUE
)
ggsave(file=paste0(projectName,"/figures/","VIDGER-MA.png"), vmaplot, width = 18, dpi = 300, units = "cm", device='png')


vplot <- vsVolcano(
    x = condControl, y = condTreated, 
    data = dds, d.factor = varInt, type = 'deseq', 
    padj = 0.05, x.lim = NULL, lfc = 1, title = TRUE, 
    legend = TRUE, grid = TRUE, data.return = FALSE
)

ggsave(file=paste0(projectName,"/figures/","VIDGER-VOLCANO.png"), vplot, width = 18, dpi = 300, units = "cm", device='png')


fpm_distribution <-vsBoxPlot(
    data = dds, d.factor = varInt, type = 'deseq', 
    title = TRUE, legend = TRUE, grid = TRUE)
ggsave(file=paste0(projectName,"/figures/","FPM_Distribution.png"), fpm_distribution, width = 18, dpi = 300, units = "cm", device='png')


scatter_plot<-vsScatterPlot(
    x = condControl, y = condTreated, 
    data = dds, type = 'deseq', d.factor = varInt, 
    title = TRUE, grid = TRUE)

ggsave(file=paste0(projectName,"/figures/","scatterplot.png"), scatter_plot, width = 18, dpi = 300, units = "cm", device='png')


vsBoxPlot<-vsBoxPlot(
    data = dds, d.factor = varInt, type = 'deseq', 
    title = TRUE, legend = TRUE, grid = TRUE
)
ggsave(file=paste0(projectName,"/figures/","boxplot.png"), vsBoxPlot, width = 18, dpi = 300, units = "cm", device='png')


vsVioDot<-vsBoxPlot(
    data = dds, d.factor = varInt, type = 'deseq', 
    title = TRUE, legend = TRUE, grid = TRUE, aes = "viosumm", fill.color = "Paired"
)
ggsave(file=paste0(projectName,"/figures/","vsVioDot.png"), vsVioDot, width = 18, dpi = 300, units = "cm", device='png')


library(limma)
library(pheatmap)
library(Biobase)
library(ggplot2)
library(gplots)
library(GEOquery)
library(reshape2)
library(plyr)
#1.get data from GEO
series <- "GSE9476"
gset <- getGEO(series,GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = "data/")
length(gset)
class(gset)
names(gset)
gset <- gset[[1]]         #platform<- "GPL96"
#if (lenght(gse)>1) {idx <- grep(platform,attr(gset,"names"))} else {idx <- 1}
#gset <- gset [[idx]]
gr <- c("CD34",rep("BM",10),rep("CD34",7),rep("AML",26),rep("PB",10),rep("CD34",10))
lenght(gr)

#2.get gene expression matrix
ex <- exprs(gset)
dim(ex)
max(ex)
min(ex)                   #ex <- log2(ex+1)
#exprs(gset) <- ex

#3.Quality control stages:
#3-1:BOXPLOT:
pdf("results/boxplot.pdf",width = 64)
boxplot(ex)
dev.off()                 #ex <- normalizequantiles(ex)
#exprs(gset) <- ex
#3-2:HEATMAP:
pdf("results/corHeatmap.pdf",width = 15, height = 15)
pheatmap(cor(ex), labels_row = gr, labels_col = gr)
dev.off()

#3-3.PRINCIPLE COMPONENT ANALYSIS:(GENES)
pc <- prcomp(ex)
pdf("results/pc.pdf")
plot(pc)
dev.off()
names(pc)
plot[pc$x[,1:2]]
dev.off()
dim(pc$x)
colnames(pc$x)
#instead:
ex.scale <- t(scale(t(ex),scale = FALSE))
pc <- prcomp(ex.scale)
pdf("results/pc_scale.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

#3-4.PRINCIPLE COMPONENT ANALYSIS:(samples)
pcr <- data.frame(pc$rottation[,1:3],Group = gr)
head(pcr)
pdf("result/pc_Ssamples.pdf")
ggplot(pcr,aes(pc1,pc2,color=Group)) + geom-point(size=3) + theme_bw()
dev.off()

#Differential expression analysis:
gr <- as.factor(gr)
gset$description <- gr
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(gr)
head(design)
fit <- lmFit(gset,design)
cont.matrix <- makeContrasts(AML-CD34, levels = design)
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2,0.01)
tT <- toptable(fit2,adjust.method = "fdr",sort.by = "B",number = Inf)
head(tT)
colnames(tT)
head(tT$Gene.ID)
tT <- subset (tT,select = c("Gene.Symbol","Gene.ID","adj.P.value","logfc"))
write.table(tT,"resulst/AML-CD34.txt",row.names = FALSE, sep = "\t",quote = FALSE)

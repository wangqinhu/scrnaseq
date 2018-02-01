library(SingleCellExperiment)
library(scater)
library(scran)


# sce config
total_num_cell <- 1920
tiss_num_cell <- 384
tiss <- c("E", "M", "A", "L", "P")
tiss_name <- c("EpiSCs", "Middle PS", "Anterior PS", "Lateral Mesoderm", "Paraxial Mesoderm")
tiss_col <- c("#00B0F0", "#00B050", "#7F7F7F", "#FE8287", "#6D1299")
cell_name <- c(paste("E", seq(1:tiss_num_cell), sep = ""),
               paste("M", seq(1:tiss_num_cell), sep = ""),
               paste("A", seq(1:tiss_num_cell), sep = ""),
               paste("L", seq(1:tiss_num_cell), sep = ""),
               paste("P", seq(1:tiss_num_cell), sep = "")
)

################
# Count loading
###############
scrnaseq_counts_file <- "../data/matrix/scRNA-seq.counts.tsv"
all.counts <- as.matrix(read.table(scrnaseq_counts_file, sep = "\t", header = TRUE, row.names = 1))
colnames(all.counts) <- cell_name

# create a sce
sce <- SingleCellExperiment(list(counts=all.counts))
dim(sce)
# set spike-in
is.spike <- grepl("^ERCC", rownames(sce))
isSpike(sce, "ERCC") <- is.spike
summary(is.spike)
# set mito
is.mito <- grepl("^mt-", rownames(sce))
summary(is.mito)


##############################
# Quality control on the cells
##############################

# calculate QC metrics
sce <- calculateQCMetrics(sce, feature_controls=list(ERCC=is.spike, Mt=is.mito))
head(colnames(colData(sce)))

# plot distributions of QC metrics
pdf("scplot/QC_metrics-v2.pdf", width = 6, height = 6)
par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))
hist(sce$total_counts/1e3, xlab="Library sizes (thousands)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce$pct_counts_ERCC, xlab="ERCC proportion (%)", 
     ylab="Number of cells", breaks=100, main="", col="grey80")
hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=100, main="", col="grey80")
dev.off()

# Removing low-quality cells based on outliers
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)
spike.drop <- isOutlier(sce$pct_counts_ERCC, nmads=3, type="higher")
sce <- sce[,!(libsize.drop | feature.drop | spike.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           BySpike=sum(spike.drop), Remaining=ncol(sce))

# PCA based on the quality metrics
pdf("scplot/PCA_QC.pdf", width = 5, height = 5)
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotPCA(sce, pca_data_input="pdata") + fontsize
dev.off()

# cell color
cell_index <- charmatch(substr(colnames(sce), 1, 1), tiss)
cell_col <- tiss_col[cell_index]


####################################
# Classification of cell cycle phase
####################################
set.seed(2018)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
assignments <- cyclone(sce, mm.pairs, gene.names=rownames(sce))

pdf("scplot/cell_cycle.pdf", width = 5, height = 5)
par(mar=c(4,4,1,1))
plot(assignments$score$G1, assignments$score$G2M, 
     xlab="G1 score", ylab="G2/M score", pch=16, cex = 0.8, col=cell_col)
legend("topright", tiss_name,
       pch=16, cex = 0.75,
       col = tiss_col,
       box.lwd = NA)
dev.off()

sce$phases <- assignments$phases
table(sce$phases)

# cell cycle of five pops

# pop1: E
sce_E <- sce[, cell_index == 1]
assignments_E <- cyclone(sce_E, mm.pairs, gene.names=rownames(sce_E))
pdf("scplot/cell_cycle_E.pdf", width = 5, height = 5)
par(mar=c(4,4,1,1))
plot(assignments_E$score$G1, assignments_E$score$G2M, 
     xlab="G1 score", ylab="G2/M score", pch=16, cex = 0.8)
dev.off()
sce_E$phases <- assignments_E$phases
table(sce_E$phases)
#

# pop2: M
sce_M <- sce[, cell_index == 2]
assignments_M <- cyclone(sce_M, mm.pairs, gene.names=rownames(sce_M))
pdf("scplot/cell_cycle_M.pdf", width = 5, height = 5)
par(mar=c(4,4,1,1))
plot(assignments_M$score$G1, assignments_M$score$G2M, 
     xlab="G1 score", ylab="G2/M score", pch=16, cex = 0.8)
dev.off()
sce_M$phases <- assignments_M$phases
table(sce_M$phases)
#

# pop3: A
sce_A <- sce[, cell_index == 3]
assignments_A <- cyclone(sce_A, mm.pairs, gene.names=rownames(sce_A))
pdf("scplot/cell_cycle_A.pdf", width = 5, height = 5)
par(mar=c(4,4,1,1))
plot(assignments_A$score$G1, assignments_A$score$G2M, 
     xlab="G1 score", ylab="G2/M score", pch=16, cex = 0.8)
dev.off()
sce_A$phases <- assignments_A$phases
table(sce_A$phases)
#

# pop4: L
sce_L <- sce[, cell_index == 4]
assignments_L <- cyclone(sce_L, mm.pairs, gene.names=rownames(sce_L))
pdf("scplot/cell_cycle_L.pdf", width = 5, height = 5)
par(mar=c(4,4,1,1))
plot(assignments_L$score$G1, assignments_L$score$G2M, 
     xlab="G1 score", ylab="G2/M score", pch=16, cex = 0.8)
dev.off()
sce_L$phases <- assignments_L$phases
table(sce_L$phases)
#

# pop5: P
sce_P <- sce[, cell_index == 5]
assignments_P <- cyclone(sce_P, mm.pairs, gene.names=rownames(sce_P))
pdf("scplot/cell_cycle_P.pdf", width = 5, height = 5)
par(mar=c(4,4,1,1))
plot(assignments_P$score$G1, assignments_P$score$G2M, 
     xlab="G1 score", ylab="G2/M score", pch=16, cex = 0.8)
dev.off()
sce_P$phases <- assignments_P$phases
table(sce_P$phases)
#

#########################################
# Examining gene-level expression metrics
#########################################

# Inspecting the most highly expressed genes
pdf("scplot/Top50gene.pdf", width = 6, height = 8)
par(mar=c(4,4,1,1))
plotQC(sce, type = "highest-expression", n=50) + fontsize
dev.off()

# Filtering out low-abundance genes

# - gene filter 1: filter out gene which average expr less than 0.01 
ave.counts <- calcAverage(sce)
pdf("scplot/gene.average_count.pdf", width = 4, height = 4)
par(mar=c(4,4,1,1))
hist(log10(ave.counts), breaks=100, main="", col="grey80", 
     xlab=expression(Log[10]~"average count"))
dev.off()

gene.keep <- ave.counts >= 0.01
filtered.sce <- sce[gene.keep,]
summary(gene.keep)

# - gene filter 2: filter out gene which expr in less than 1% of the total cell
filtered.ave.counts <- calcAverage(filtered.sce)
num.cells <- nexprs(filtered.sce, byrow=TRUE)
pdf("scplot/gene_num.percell.pdf", width = 4, height = 4)
par(mar=c(4,4,1,1))
smoothScatter(log10(filtered.ave.counts), num.cells, ylab="Number of cells", 
              xlab=expression(Log[10]~"average count"))
dev.off()

to.keep <- num.cells > total_num_cell/100
sce <- filtered.sce[to.keep,]
summary(to.keep)

#######################################
# Normalization of cell-specific biases
#######################################
sce <- computeSumFactors(sce)
summary(sizeFactors(sce))

pdf("scplot/sizefactor_libsize.pdf", width = 4, height = 4)
par(mar=c(4,4,1,1))
plot(sizeFactors(sce), sce$total_counts/1e6, log="xy",
     ylab="Library size (millions)", xlab="Size factor")
dev.off()

sce <- computeSpikeFactors(sce, type="ERCC", general.use=FALSE)

sce <- normalize(sce)

############################################
# Checking for confounding technical factors
############################################
pdf("scplot/var_counts_ercc.pdf", width = 6, height = 6)
par(mar=c(4,4,1,1))
plotExplanatoryVariables(sce, variables=c("total_counts_ERCC", 
                                          "log10_total_counts_ERCC")) + fontsize
dev.off()

##################################################
# Modelling the technical noise in gene expression
##################################################
var.fit.nospike <- trendVar(sce, parametric=TRUE, use.spikes=FALSE, span=0.2)
var.out.nospike <- decomposeVar(sce, var.fit.nospike)

pdf("scplot/trendfit.pdf", width = 4, height = 4)
par(mar=c(4,4,1,1))
plot(var.out.nospike$mean, var.out.nospike$total, pch=16, cex=0.6, 
     xlab="Mean log-expression", ylab="Variance of log-expression")
curve(var.fit.nospike$trend(x), col="dodgerblue", lwd=2, add=TRUE)
cur.spike <- isSpike(sce)
points(var.out.nospike$mean[cur.spike], var.out.nospike$total[cur.spike], col="red", pch=16)
dev.off()

chosen.genes <- order(var.out.nospike$bio, decreasing=TRUE)[1:10]
pdf("scplot/top10_vioplot.pdf", width = 6, height = 6)
par(mar=c(4,4,1,1))
plotExpression(sce, features=rownames(var.out.nospike)[chosen.genes]) + fontsize
dev.off()

# Highly variable genes (HVGs)
top.hvgs <- order(var.out.nospike$bio, decreasing=TRUE)
head(var.out.nospike[top.hvgs,], 100)

# Denoising expression values using PCA
sce <- denoisePCA(sce, technical=var.fit.nospike$trend) 
dim(reducedDim(sce, "PCA"))

# PCA plot
pdf("scplot/feature.PCA.pdf", width = 6, height = 6)
par(mar=c(4,4,1,1))
plotReducedDim(sce, use_dimred="PCA", ncomponents=3, colour_by="total_features") + fontsize
dev.off()

# tSNE plot
pdf("scplot/feature.tSNE.pdf", width = 16, height = 5)
out5 <- plotTSNE(sce, use_dimred="PCA", perplexity=5, colour_by="total_features", 
                 rand_seed=100) + fontsize + ggtitle("Perplexity = 5")
out10 <- plotTSNE(sce, use_dimred="PCA", perplexity=10, colour_by="total_features",
                  rand_seed=100) + fontsize + ggtitle("Perplexity = 10")
out20 <- plotTSNE(sce, use_dimred="PCA", perplexity=20, colour_by="total_features",
                  rand_seed=100) + fontsize + ggtitle("Perplexity = 20")
multiplot(out5, out10, out20, cols=3)
dev.off()

sessionInfo()

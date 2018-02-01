library(gplots)
library(SingleCellExperiment)
library(scater)
library(scran)
library(Rtsne)


load("sce.RData")

cell_mark <- c(8,15:18)

#x <- logcounts(sce)
x <- counts(sce)

# top 10% hvg
keep <- rownames(var.out.nospike[top.hvgs,][1:as.integer(dim(sce)[1]/10),])
y <- x[keep,]

sc_unique <- unique(t(y)) # Remove duplicates
sc_matrix <- as.matrix(sc_unique)
set.seed(2017) # Set a seed if you want reproducible results
tsne_out1 <- Rtsne(sc_matrix) # Run TSNE

tiss_name<-c("EpiSCs", "Middle PS", "Anterior PS", "Lateral Mesoderm", "Paraxial Mesoderm")

# Show the objects in the 2D tsne representation
pdf("scplot/AAAA-2D-tSNE1.pdf", 5, 5)
par(mar=c(4,4,1,1))
plot(tsne_out1$Y,col=cell_col, pch=16, cex = 0.75,
     xlab = "t-SNE 1", ylab= "t-SNE 2",
     xlim = c(-40,40),
     ylim = c(-40,40)
     )
legend("topleft", tiss_name,
       pch=16, cex = 0.75,
       pch = tiss_col,
       box.lwd = NA)
dev.off()

# Show the objects in the 2D tsne representation, marker gene color map
makr_gene <- "ENSMUSG00000062327"
mark_gene_counts <- x[makr_gene,]
pdf("scplot/AAAA-2D-tSNE1-Pbx1.pdf", 5, 5)
par(mar=c(4,4,1,1))
plot(tsne_out1$Y,
     col=greenred(max(mark_gene_counts))[mark_gene_counts],
     pch=cell_mark[cell_index], cex = 0.75,
     xlab = "t-SNE 1", ylab= "t-SNE 2",
     xlim = c(-40,40),
     ylim = c(-40,40)
)
legend("topleft", tiss_name,
       pch=cell_mark, cex = 0.75,
       box.lwd = NA)
dev.off()

# Using a dist object
tsne_out2 <- Rtsne(dist(sc_matrix))

pdf("scplot/AAAA-2D-tSNE2.pdf", 5, 5)
par(mar=c(4,4,1,1))
plot(tsne_out2$Y,col=cell_col, pch=16, cex = 0.75,
     xlab = "t-SNE 1", ylab= "t-SNE 2",
     xlim = c(-40,40),
     ylim = c(-40,40)
)
legend("topleft", tiss_name,
       pch=16, cex = 0.75,
       col = tiss_col,
       box.lwd = NA)
dev.off()


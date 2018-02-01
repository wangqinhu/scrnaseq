library(gplots)
library(SingleCellExperiment)
library(scater)
library(scran)
library(Rtsne)
library(scatterplot3d)
library(rgl)

load("sce.RData")

x <- logcounts(sce)

sc_unique <- unique(t(x)) # Remove duplicates
sc_matrix <- as.matrix(sc_unique)
set.seed(2017) # Set a seed if you want reproducible results
tsne_out <- Rtsne(dist(sc_matrix), dims = 3) # Run TSNE

pdf("scplot/3D-tSNE_distance.pdf", width = 6, height = 6)
scatterplot3d(x = tsne_out$Y[,1], y = tsne_out$Y[,2], z = tsne_out$Y[,3],
              pch = 16, cex.symbols = 0.5,
              xlab = "t-SNE 1", ylab= "t-SNE 2", zlab = "t-SNE 3",
              color = cell_col)
legend("topleft", tiss_name,
       pch=16, cex = 0.75,
       col = tiss_col,
       bg = NULL,
       box.lwd = NA)
dev.off()

# free rotate
plot3d(x=tsne_out$Y[,1], y=tsne_out$Y[,2], z=tsne_out$Y[,3],
       pch = 16,
       xlab = "t-SNE 1", ylab= "t-SNE 2", zlab = "t-SNE 3",
       col = cell_col)


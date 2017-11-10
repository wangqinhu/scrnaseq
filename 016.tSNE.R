library(Rtsne)

tiss_col<-c("#00B0F0", "#00B050", "#7F7F7F", "#FE8287", "#6D1299")
cell_num<-c(327, 335, 361, 330, 302)

cell_col<-c(rep(tiss_col[1], cell_num[1]),
            rep(tiss_col[2], cell_num[2]),
            rep(tiss_col[3], cell_num[3]),
            rep(tiss_col[4], cell_num[4]),
            rep(tiss_col[5], cell_num[5])
)

x<-read.delim(file="tSNE.matrix.tsv", row.names=1, stringsAsFactors=FALSE, as.is=TRUE, sep ="\t")


sc_unique <- unique(x)
sc_matrix <- as.matrix(sc_unique)
set.seed(2017)
tsne_out <- Rtsne(sc_matrix)

tiss_name<-c("EpiSCs", "Middle PS", "Anterior PS", "Lateral Mesoderm", "Paraxial Mesoderm")

pdf("t-SNE.pdf", 5, 5)
par(mar=c(4,4,1,1))
plot(tsne_out$Y,col=cell_col, pch=16, cex = 0.75,
     xlab = "t-SNE 1", ylab= "t-SNE 2",
     xlim = c(-40,40),
     ylim = c(-40,40)
     )
legend("topleft", tiss_name,
       pch=16, cex = 0.75,
       col = tiss_col,
       box.lwd = NA)
dev.off()

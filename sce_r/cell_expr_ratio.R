data <- read.table("../data/matrix/scRNA-seq.counts.tsv", header=TRUE, row.names=1)

x <- as.matrix(data)

y<-rep(0,100)
for (i in 1:100) {
  keep <- rowSums(x>=1) >= dim(data)[2]/5 * (i/100)
  y[i] <- dim(x[keep,])
}
pdf("scplot/cell_expr_ratio.pdf", 4, 4)
par(mar=c(4,4,1,1))
plot(1:100,y, pch=16, col=heat.colors(100), xlab = "Percentage (%)", ylab="Number of genes", type = "p")
dev.off()
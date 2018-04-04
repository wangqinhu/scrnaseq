#options(scipen=5)

rm(list=ls())
x<-read.table("data/matrix/guo2.stat.tsv", header = T, row.names = 1)

pdf("scRNAseq_stats.guo2.pdf", width = 10, height = 3)
layout(matrix(c(1:6),1,6,byrow = TRUE))
par(mar=c(2.5,4,2,0.2))

x$tissue<-ordered(x$tissue, levels=c("B", "WT", "15", "28"))
tiss_col<-c("red", "#1E90FF", "#00BFFF", "#97FFFF")
#tiss_name<-c("EpiSCs", "Middle PS", "Anterior PS", "Lateral Mesoderm", "Paraxial Mesoderm")

boxplot(x$barcoded~x$tissue, col=tiss_col, outpch=16, outcol="grey", frame = FALSE, ylim=c(0,800000), main="Barcoded", ylab="Number of reads")
boxplot(x$short~x$tissue, col=tiss_col, outpch=16, outcol="grey", frame = FALSE, ylim=c(0,800000), main="Length < 30 nt")
boxplot(x$clean~x$tissue, col=tiss_col, outpch=16, outcol="grey", frame = FALSE, ylim=c(0,800000), main="Clean")
boxplot(x$mapped~x$tissue, col=tiss_col, outpch=16, outcol="grey", frame = FALSE, ylim=c(0,800000), main="Mapped")
boxplot(x$genic_unique~x$tissue, col=tiss_col, outpch=16, outcol="grey", frame = FALSE, ylim=c(0,80000), main="Genic-unique")
boxplot(x$gene~x$tissue, col=tiss_col, outpch=16, outcol="grey", frame = FALSE, ylim=c(0,8000), main="Number of genes")

#legend("topright", tiss_name, pch=16, cex=0.9, col = tiss_col, box.lwd = NA)
dev.off()

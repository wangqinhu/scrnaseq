#options(scipen=5)

rm(list=ls())
x<-read.table("../data/matrix/scRNAseq_stats.txt", header = T, row.names = 1)

pdf("scplot/scRNAseq_stats.pdf", width = 10, height = 3)
layout(matrix(c(1:6),1,6,byrow = TRUE))
par(mar=c(2.5,4,2,0.2))

x$tissue<-ordered(x$tissue, levels=c("E", "M", "A", "L", "P"))

tiss_col<-c("#00B0F0", "#00B050", "#7F7F7F", "#FE8287", "#6D1299")
#tiss_name<-c("EpiSCs", "Middle PS", "Anterior PS", "Lateral Mesoderm", "Paraxial Mesoderm")

boxplot(x$barcoded~x$tissue, col=tiss_col, outpch=16, outcol="grey", frame = FALSE, ylim=c(0,800000), main="Barcoded", ylab="Number of reads")
boxplot(x$short~x$tissue, col=tiss_col, outpch=16, outcol="grey", frame = FALSE, ylim=c(0,800000), main="Length < 30 nt")
boxplot(x$clean~x$tissue, col=tiss_col, outpch=16, outcol="grey", frame = FALSE, ylim=c(0,800000), main="Clean")
boxplot(x$mapped~x$tissue, col=tiss_col, outpch=16, outcol="grey", frame = FALSE, ylim=c(0,800000), main="Mapped")
boxplot(x$genic_unique~x$tissue, col=tiss_col, outpch=16, outcol="grey", frame = FALSE, ylim=c(0,80000), main="Genic-unique")
boxplot(x$gene~x$tissue, col=tiss_col, outpch=16, outcol="grey", frame = FALSE, ylim=c(0,8000), main="Number of genes")

#legend("topright", tiss_name, pch=16, cex=0.9, col = tiss_col, box.lwd = NA)
dev.off()

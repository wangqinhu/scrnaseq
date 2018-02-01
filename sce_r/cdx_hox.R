library(gplots)
library(SingleCellExperiment)
library(scater)
library(scran)

load("sce.RData")
load("hc.RData")

x <- logcounts(sce)

cdx<-read.table("cdx.txt")
hox<-read.table("hox.txt")

cdx_index <- charmatch(as.vector(cdx$V1), rownames(x))
hox_index <- charmatch(as.vector(hox$V1), rownames(x))
cdx_mat <- x[cdx_index,]
hox_mat <- x[hox_index,]
rownames(cdx_mat)<-as.vector(cdx$V2)
rownames(hox_mat)<-as.vector(hox$V2)
filtered.cdx <- na.omit(as.matrix(cdx_mat))
filtered.hox <- na.omit(as.matrix(hox_mat))

# cdx heatmap
y <- filtered.cdx

hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="ward.D2")

myro <- cutree(hr, h=max(hr$height)/2, 10);
myrohr <- rainbow(length(unique(myro)))
myrohr <- myrohr[as.vector(myro)]

cyre=colorRampPalette(c("cyan","white","red"),space="rgb")
myheatcol <- cyre(4)

pdf("scplot/cdx_heatmap.pdf", width = 5, height = 5)

z<-heatmap.2(y, 
             Rowv=as.dendrogram(hr),
             Colv=as.dendrogram(hc), 
             col=myheatcol,
             scale = "row",
             density.info="none",
             trace="none",
             labCol = NA,
             margins=c(0.6,6),
             ColSideColors=cell_col,
             RowSideColors=myrohr,
             key.title = NA,
             key.xlab = NA,
             key.xtickfun=function() {
               cex <- par("cex")*par("cex.axis")
               side <- 1
               line <- -1.5
               col <- par("col.axis")
               font <- par("font.axis")
               mtext("Low", side=side, at=-0.05, adj=1,
                     line=line, cex=cex, col=col, font=font)
               mtext("High", side=side, at=1.05, adj=0,
                     line=line, cex=cex, col=col, font=font)
               return(list(labels=FALSE, tick=FALSE))
             },
             key.par=list(
               par(mar=c(0.1,2,8.9,1)), cex=0.6
             )
)
par(xpd=T)
legend(-0.2,1.25, tiss_name,
       pch=16, cex = 0.75,
       col = tiss_col,
       box.lwd = NA)
dev.off()

# hox heatmap
y <- filtered.hox

hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="ward.D2")

myro <- cutree(hr, h=max(hr$height)/2, 10);
myrohr <- rainbow(length(unique(myro)))
myrohr <- myrohr[as.vector(myro)]


cyre=colorRampPalette(c("cyan","white","red"),space="rgb")
myheatcol <- cyre(4)

pdf("scplot/hox_heatmap.pdf", width = 5, height = 5)

z<-heatmap.2(y, 
             Rowv=as.dendrogram(hr),
             Colv=as.dendrogram(hc), 
             col=myheatcol,
             scale = "row",
             density.info="none",
             trace="none",
             labCol = NA,
             margins=c(0.6,6),
             ColSideColors=cell_col,
             RowSideColors=myrohr,
             key.title = NA,
             key.xlab = NA,
             key.xtickfun=function() {
               cex <- par("cex")*par("cex.axis")
               side <- 1
               line <- -1.5
               col <- par("col.axis")
               font <- par("font.axis")
               mtext("Low", side=side, at=-0.05, adj=1,
                     line=line, cex=cex, col=col, font=font)
               mtext("High", side=side, at=1.05, adj=0,
                     line=line, cex=cex, col=col, font=font)
               return(list(labels=FALSE, tick=FALSE))
             },
             key.par=list(
               par(mar=c(0.1,2,8.9,1)), cex=0.6
             )
)
par(xpd=T)
legend(-0.2,1.25, tiss_name,
       pch=16, cex = 0.75,
       col = tiss_col,
       box.lwd = NA)
dev.off()
library(gplots)
library(SingleCellExperiment)
library(scater)
library(scran)

load("sce.RData")


x <- logcounts(sce)
# top 10% hvg
keep <- rownames(var.out.nospike[top.hvgs,][1:as.integer(dim(sce)[1]/10),])
y <- x[keep,]

hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="ward.D2")
hc <- hclust(as.dist(1-cor(y, method="spearman")), method="ward.D2") 


myco <- cutree(hc, h=max(hc$height)/2, 10);
mycohc <- rainbow(length(unique(myco)))
mycohc <- mycohc[as.vector(myco)]

myro <- cutree(hr, h=max(hr$height)/2, 10);
myrohr <- rainbow(length(unique(myro)))
myrohr <- myrohr[as.vector(myro)]


cyre=colorRampPalette(c("cyan","white","red"),space="rgb")
myheatcol <- cyre(4)

pdf("scplot/gene_heatmap.pdf", width = 5, height = 5)

z<-heatmap.2(y, 
          Rowv=as.dendrogram(hr),
          Colv=as.dendrogram(hc), 
          col=myheatcol,
          scale = "row",
          density.info="none",
          trace="none",
          labRow = NA,
          labCol = NA,
          margins=c(0.6,0.6),
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
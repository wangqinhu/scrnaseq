library(gplots)
library(SingleCellExperiment)
library(scater)
library(scran)

load("sce.RData")


x <- logcounts(sce)
# top 10% hvg
keep <- rownames(var.out.nospike[top.hvgs,][1:as.integer(dim(sce)[1]/10),])
y <- x[keep,]

corr<-cor(y, method = c("spearman"))
cyre=colorRampPalette(c("cyan","white","red"),space="rgb")

par(mar=c(1,1,1,1))
pdf("scplot/cell_heatmap.pdf", 5, 5)
z<-heatmap.2(corr, RowSideColors = cell_col, ColSideColors = cell_col,
          trace="none", density.info="none",
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
          ),
          labRow = NA, labCol = NA,
          margins=c(0.6,0.6),
          col=cyre(4), scale="row")
par(xpd=T)
legend(-0.2,1.25, tiss_name,
       pch=16, cex = 0.75,
       col = tiss_col,
       box.lwd = NA)
dev.off()

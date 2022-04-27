## ----fig.height=5, fig.width=5, message=FALSE, warning=FALSE------------------
library(DESeq2)
library(lfc)
data(airway, package="airway")
A <- assay(airway)[,2]
B <- assay(airway)[,1]

ll <- PsiLFC(A,B)
plot(ecdf(ll),xlim=c(-1,1),xlab="Log2 fold change treated/untreated",ylab="Cumulative frequency",main="Cell N61311",col='blue')
lines(ecdf(CenterMedian(log2((A+1)/(B+1)))),col='red')

## ----fig.width=5, fig.height=5------------------------------------------------
plot(ecdf(ll[A>0 | B>0]),xlim=c(-1,1),xlab="Log2 fold change treated/untreated",ylab="Cumulative frequency",main="Cell N61311",col='blue')
lines(ecdf(CenterMedian(log2((A+1)/(B+1))[A>0 | B>0])),col='red')

## -----------------------------------------------------------------------------
head(PsiLFC.se(airway,contrast=c("dex","untrt","trt")))

## ----message=FALSE, warning=FALSE---------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = assay(airway),colData = colData(airway),design= ~ dex)
dds <- DESeq(dds)
res <- results(dds, contrast=c("dex","untrt","trt"),cre=TRUE)
head(res)


## ----fig.height=5, fig.width=5, message=FALSE, warning=FALSE------------------
library(DESeq2)
library(lfc)

# it is overly complicated to read a gzipped table from the internet. This is how it is done:
t <- read.delim(textConnection(readLines(gzcon(url("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE160764&format=file&file=GSE160764%5FRNA%5FMEF%5FOla%5FIFNb%5FCounttable%5FRaw%2Etxt%2Egz")))))

A <- t[,4]
B <- t[,2]

ll <- PsiLFC(A,B)
plot(ecdf(ll),xlim=c(-1,1),xlab="Log2 fold change treated/untreated",ylab="Cumulative frequency",col='blue')
lines(ecdf(CenterMedian(log2((A+1)/(B+1)))),col='red')

## ----fig.width=5, fig.height=5------------------------------------------------
plot(ecdf(ll[A>0 | B>0]),xlim=c(-1,1),xlab="Log2 fold change treated/untreated",ylab="Cumulative frequency",col='blue')
lines(ecdf(CenterMedian(log2((A+1)/(B+1))[A>0 | B>0])),col='red')

## ----message=FALSE, warning=FALSE---------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = t[,2:5],colData = data.frame(IFN=c("0h","0h","1h","1h")),design= ~ IFN)
dds <- DESeq(dds)
res <- results(dds, contrast=c("IFN","1h","0h"),cre=TRUE)
head(res)


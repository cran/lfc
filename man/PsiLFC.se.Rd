% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psi_lfc.R
\name{PsiLFC.se}
\alias{PsiLFC.se}
\title{Psi LFC effect size estimator}
\usage{
PsiLFC.se(se, contrast, cre = FALSE)
}
\arguments{
\item{se}{SummarizedExperiment object}

\item{contrast}{Vector of length 3 (<name>,<A>,<B>)}

\item{cre}{Compute credible intervals as well? (can also be a vector of quantiles)}
}
\value{
Either a vector containing the estimates, or a data frame containing
    the credible interval as well
}
\description{
Computes the optimal effect size estimate and credible intervals if needed
for a Bioconductor SummarizedExperiment object
}
\examples{
\dontrun{
   data(airway, package="airway")
   head(PsiLFC.se(airway,contrast=c("dex","untrt","trt")))
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rem.R
\name{rem}
\alias{rem}
\title{Exclude Genes with Low Expression}
\usage{
rem(x, missing.perc, const.int)
}
\arguments{
\item{x}{The expression matrix or data frame with genes as rows and samples as columns.}

\item{missing.perc}{A cutoff to remove genes with zero expression across samples.}

\item{const.int}{A pseudocount is added to the TPM (Transcripts Per Million) values before performing a log transformation.}
}
\value{
Processed expression data after removing low-expressed genes.
}
\description{
Exclude genes with low or zero expression.
}

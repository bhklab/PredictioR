% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneSurvCont.R
\name{geneSurvCont}
\alias{geneSurvCont}
\title{Fit Proportional Hazards Regression Model for Genes: Continuous Expression Variable}
\usage{
geneSurvCont(
  dat.icb,
  clin = NULL,
  time.censor,
  missing.perc,
  const.int = 0.001,
  n.cutoff,
  feature,
  study,
  surv.outcome,
  cancer.type,
  treatment
)
}
\arguments{
\item{dat.icb}{A MultiAssayExperiment (MAE) object, SummarizedExperiment (SE) object, or a data frame or matrix of gene expression data.}

\item{clin}{If dat.icb is a data frame or matrix, then it contains clinical data (as data frame or matrix). By default, it is NULL.}

\item{time.censor}{Possible censoring in months.}

\item{missing.perc}{A cutoff to remove genes with zero expression across samples.}

\item{const.int}{A pseudocount is added to the TPM (Transcripts Per Million) values before performing a log transformation.}

\item{n.cutoff}{Minimum number of samples included in the association analysis.}

\item{feature}{A vector of character strings for selected features.}

\item{study}{Name of study.}

\item{surv.outcome}{Overall survival (OS) or progression-free survival (PFS).}

\item{cancer.type}{Name of the cancer type for the given study.}

\item{treatment}{Name of the treatment for the given study.}
}
\value{
A subset of results using an object of class 'coxph' representing the fit.
Outcome: Immunotherapy time-to-event outcomes including overall survival (OS) and progression-free survival (PFS).
Gene: Name of selected genes.
Study: Name of study.
Coef: Estimate of treatment effect i.e., log hazard ratio.
SE: Standard error of treatment estimate.
N: Number of samples.
Pval: Estimated p-value.
Cancer_type: A character shows the cancer type.
Treatment: A character shows the treatment type.
}
\description{
Fits a Cox proportional hazards regression model with continuous expression data for each gene using the counting process formulation of Andersen and Gill.
}
\examples{
Assess the association between selected genes and OS in immunotherapy.
geneSurvCont(dat.icb = ICB_small_Mariathasan,
             time.censor = 36,
             missing.perc = 0.5,
             const.int = 0.001,
             n.cutoff = 15,
             feature = c('CXCL9', 'CXCL10', 'TIGIT', 
                         'CD83', 'STAT1', 'CXCL11',
                         'CXCL13', 'CD8A', 'CTLA4'),
             study = 'ICB_Mariathasan',
             surv.outcome = 'OS',
             cancer.type = 'Bladder',
             treatment = 'PD-1/PD-L1')
             
}

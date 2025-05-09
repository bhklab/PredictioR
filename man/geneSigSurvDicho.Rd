% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneSigSurvDicho.R
\name{geneSigSurvDicho}
\alias{geneSigSurvDicho}
\title{Fit Proportional Hazards Regression Model for Gene Signature: Dichotomous Variable}
\usage{
geneSigSurvDicho(
  dat.icb,
  clin = NULL,
  geneSig,
  time.censor,
  n.cutoff,
  n0.cutoff,
  n1.cutoff,
  study,
  surv.outcome,
  sig.name,
  method = "median",
  var.type,
  cancer.type,
  treatment
)
}
\arguments{
\item{dat.icb}{A MultiAssayExperiment (MAE) object, SummarizedExperiment (SE) object, or a data frame or matrix of gene expression data.}

\item{clin}{If dat.icb is a data frame or matrix, then it contains clinical data (as data frame or matrix). By default, it is NULL.}

\item{geneSig}{A numeric vector of computed signature score.}

\item{time.censor}{Possible censoring in months.}

\item{n.cutoff}{Minimum number of samples included in the association analysis.}

\item{n0.cutoff}{Minimum number of samples with status 0.}

\item{n1.cutoff}{Minimum number of samples with status 1.}

\item{study}{Name of study.}

\item{surv.outcome}{Overall survival (OS) or progression-free survival (PFS).}

\item{sig.name}{Name of signature.}

\item{method}{The default method to convert a continuous variable into a dichotomous variable is the 'median' method. The first quartile (Q1) and third quartile (Q3) can also be applied.}

\item{var.type}{If the variable is dichotomous (by default), then var.type is TRUE.}

\item{cancer.type}{Name of the cancer type for the given study.}

\item{treatment}{Name of the treatment for the given study.}
}
\value{
A subset of results using an object of class 'coxph' representing the fit.
Outcome: Immunotherapy time-to-event outcomes including overall survival (OS) and progression-free survival (PFS).
Gene: Name of selected signature.
Study: Name of study.
Coef: Estimate of treatment effect i.e., log hazard ratio.
SE: Standard error of treatment estimate.
N: Number of samples.
Pval: Estimated p-value.
Cancer_type: A character shows the cancer type.
Treatment: A character shows the treatment type.
}
\description{
Fits a Cox proportional hazards regression model with dichotomous signature data using the counting process formulation of Andersen and Gill.
}
\examples{
Assess the association of EMT Thompson signature and OS. 
sig <- geneSigMean(dat.icb = ICB_small_Mariathasan, 
                   sig = EMT_Thompson,
                   sig.name = 'EMT_Thompson',
                   missing.perc = 0.5,
                   const.int = 0.001,
                   n.cutoff = 15,
                   sig.perc = 0.8, 
                   study = 'ICB_Mariathasan')
            
geneSigSurvDicho(dat.icb = ICB_small_Mariathasan,
                 geneSig = sig,
                 time.censor = 36,
                 n.cutoff = 15,
                 n0.cutoff = 5,
                 n1.cutoff = 5,
                 study =  'ICB_Mariathasan',
                 surv.outcome = 'OS',
                 sig.name = 'EMT_Thompson',
                 method = 'median',
                 var.type = FALSE,
                 cancer.type = 'Bladder',
                 treatment = 'PD-1/PD-L1')
                
}

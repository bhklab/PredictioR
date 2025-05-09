% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneSigLogReg.R
\name{geneSigLogReg}
\alias{geneSigLogReg}
\title{Fit Logistic Regression Model for Gene Signature: Continuous Variable}
\usage{
geneSigLogReg(
  dat.icb,
  clin = NULL,
  geneSig,
  n.cutoff,
  study,
  sig.name,
  n0.cutoff,
  n1.cutoff,
  cancer.type,
  treatment
)
}
\arguments{
\item{dat.icb}{A MultiAssayExperiment (MAE) object, SummarizedExperiment (SE) object, or a data frame or matrix of gene expression data.}

\item{clin}{If dat.icb is a data frame or matrix, then it contains clinical data (as data frame or matrix). By default, it is NULL.}

\item{geneSig}{A numeric vector of computed signature score.}

\item{n.cutoff}{Minimum number of samples included in the association analysis.}

\item{study}{Name of study.}

\item{sig.name}{Name of signature.}

\item{n0.cutoff}{Minimum number of non-responders in the analysis.}

\item{n1.cutoff}{Minimum number of responders in the analysis.}

\item{cancer.type}{Name of the cancer type for the given study.}

\item{treatment}{Name of the treatment for the given study.}
}
\value{
A subset of results using an object of class logistic regression representing the fit.
Outcome: Immunotherapy response outcome i.e., R (responder) and NR (non-responder).
Gene: Name of selected signature.
Study: Name of study.
Coef: Estimate of treatment effect i.e., log odds ratio.
SE: Standard error of treatment estimate.
N: Number of samples.
Pval: Estimated p-value.
Cancer_type: A character shows the cancer type.
Treatment: A character shows the treatment type.
}
\description{
Fits a logistic regression model for continuous signature data.
}
\examples{
Assess the association of M1 Hwang signature and response. 
sig <- geneSigssGSEA(dat.icb = ICB_small_Mariathasan, 
                     sig = M1_Hwang,
                     sig.name = 'M1_Hwang',
                     missing.perc = 0.5,
                     const.int = 0.001,
                     n.cutoff = 15,
                     sig.perc = 0.8, 
                     study = 'ICB_Mariathasan')
            
geneSigLogReg(dat.icb = ICB_small_Mariathasan,
              geneSig = sig,
              n.cutoff = 15,
              study =  'ICB_Mariathasan',
              sig.name = 'M1_Hwang',
              n0.cutoff = 10,
              n1.cutoff = 10,
              cancer.type = 'Bladder',
              treatment = 'PD-1/PD-L1')
                
}

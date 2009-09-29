\name{prostate}
\encoding{latin1}
\alias{prostate}

\docType{data}

\title{Prostate cancer data, a subset data set from Stephenson et al. (2005) study}

\description{Clinical and a subset of gene expression data from the Stephenson \emph{et al.} (2005) study.}

\usage{data(prostate)}

\format{
A list containing the following components:
  \describe{
	\item{\code{type}}{vector of lenth 79 indicating the class of the patients (0 = recurrence, 1 = no recurrence).}
	\item{\code{cont}}{data matrix with 79 rows and 500 columns. 
		The gene expression of 500 randomly sampled transcripts (for memory allocation reasons, see details).}
	\item{\code{indep}}{data matrix with 79 rows and 8 columns. The measurements of 8 clinical variables. The discrete data are suitable for a \code{MElogreg} or \code{MEindep} model in the mixture of experts methodology.}
	\item{\code{loc}}{data matrix with 79 rows and 7 columns. The measurements of 7 clinical variables. The discrete data are suitable for a \code{MEloc} model in the mixture of experts methodology.}
	\item{\code{loc.ind}}{indicates the location variable.}

  }
}

\details{
The data set from Stephenson \emph{et al.} (2005) was built from tissue samples obtained from 79 patients all treated by radical prostatectomy. There were 37  samples which were classified as recurrent and 42 as non-recurrent primary prostate tumor. Samples were snap frozen and gene expression analysis was carried out using the Affymetrix U133A human gene array which has 22,283 features. After a prefiltering step, the analyzed data set contained 7,884 features. The clinical data and microarray data were measured on the same set of 79 patients.

  \itemize{
    \item clinical data: include 8 clinical variables
    \item microarray data: measure the expression of a subset of 500 randomly sampled transcripts.
  }

For the location model, variables 'semi-vesicle invasion' and 'lymph node involvement' were merged into a single categorical variable (called the location variable).

}

\source{
The data set was obtained upon request to the authors of the study. The original data with 7,884 transcripts can be downloaded as an .RData file from \url{http://www.math.univ-toulouse.fr/~lecao/package.html} 
}

\references{
Stephenson, A.J., Smith, A., Kattan, M.W., Satagopan, J., Reuter, V.E., Scardino, P.T. and Gerald, W.L. (2005). Integration of gene expression profiling and clinical variables to predict prostate carcinoma recurrence after radical prostatectomy. \emph{Cancer}, \bold{104}, 2, 290-298. 


Hunt, L. and Jorgensen, M. (1999). Mixture model clustering using the MULTIMIX program. \emph{Australian \& New Zealand Journal of Statistics}, \bold{41}, 2, 154--171. 
}

\keyword{datasets}

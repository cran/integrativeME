\name{breast}
\encoding{latin1}
\alias{breast}

\docType{data}

\title{Breast cancer data, a subset data set from van de Vijver et al. (2002) study}
\description{Clinical and a subset of gene expression data from the van de Vijver \emph{et al.} (2002) study.}

\usage{data(breast)}

\format{
A list containing the following components:
  \describe{
	\item{\code{type}}{vector of lenth 256 indicating the class of the patients (0 = recurrence, 1 = no recurrence).}
	\item{\code{cont}}{data matrix with 256  rows and 500 columns. 
		The gene expression of 500 randomly sampled transcripts (for memory allocation reasons, see details).}
	\item{\code{indep}}{data matrix with 256  rows and 8 columns. The measurements of 8 clinical variables. The discrete data are suitable for a \code{MElogreg} or \code{MEindep} model in the mixture of experts methodology.}
	\item{\code{loc}}{data matrix with 256 rows and 7 columns. The measurements of 7 clinical variables. The discrete data are suitable for a \code{MEloc} model in the mixture of experts methodology.}
	\item{\code{loc.ind}}{indicates the location variable.}

  }
}

\details{
The data set from van de Vijver \emph{et al.} (2002) contains gene expression of tumors from 256 patients who were all treated by modified radical mastectomy or breast-conserving surgery. The authors also included some patients from the Van 't Veer \emph{et al.} (2002) study and the censored patients were removed. The data were preprocessed and filtered to obtain 5,537 genes spotted on Agilent Hu25K microarrays. Eight prognostic factors were available in the clinical data and categorized as indicated by the authors.
  \itemize{
    \item clinical data: include 8 clinical variables
    \item microarray data: measure the expression of a subset of 500 randomly chosen transcripts.
  }

For the location model, variables 'posnode' and 'chemotherapy' were merged into a single categorical variable (called the location variable).

}

\source{
See website from the refered article. The original data with 5,537 transcripts can be downloaded as an .RData file from \url{http://www.math.univ-toulouse.fr/~lecao/package.html}
}

\references{
van de Vijver, M.J., He, Y.D., van' t Veer, L.J., Dai, H., Hart, A.A.M., Voskuil, D.W., Schreiber, G.J., Peterse, J.L., Roberts, C., Marton, M.J. and others (2002). A gene-expression signature as a predictor of survival in breast cancer. \emph{New England Journal of Medicine}, \bold{347}, 25, 1999--2009. 


Hunt, L. and Jorgensen, M. (1999). Mixture model clustering using the MULTIMIX program. \emph{Australian \& New Zealand Journal of Statistics}, \bold{41}, 2, 154--171. 
}

\keyword{datasets}

\name{cns}
\encoding{latin1}
\alias{cns}

\docType{data}

\title{Central Nervous System cancer data, a subset data set from Pomeroy et al. (2002) study}
\description{Clinical and a subset of gene expression data from the Pomeroy \emph{et al.} (2002) study.}

\usage{data(breast)}

\format{
A list containing the following components:
  \describe{
	\item{\code{type}}{vector of lenth 60 indicating the class of the patients (0 = dead, 1 = alive).}
	\item{\code{cont}}{data matrix with 60  rows and 500 columns. 
		The gene expression of 500 randomly sampled transcripts (for memory allocation reasons, see details).}
	\item{\code{indep}}{data matrix with 60  rows and 5 columns. The measurements of 5 clinical variables. The discrete data are suitable for a \code{MElogreg} or \code{MEindep} model in the mixture of experts methodology.}
	\item{\code{loc}}{data matrix with 60 rows and 4 columns. The measurements of 4 clinical variables. The discrete data are suitable for a \code{MEloc} model in the mixture of experts methodology.}
	\item{\code{loc.ind}}{indicates the location variable.}
	\item{\code{gene.name}}{gives more information on the genes in the \code{cont} data set.}

  }
}

\details{
Medulloblastomas are embryonal tumors of the central nervous system (CNS). Pomeroy \emph{et al.} (2002)  investigated this malignant brain tumor of childhood as the response of therapy is difficult to predict. The biopsies of 60 patients were obtained before they received any treatment. The expression level of 7,128 genes were available, as well as five clinical variables
  \itemize{
    \item clinical data: include 5 clinical variables
    \item microarray data: measure the expression of a subset of 500 randomly sampled transcripts.
  }

For the location model, variables 'Age' and 'Cx' were merged into a single categorical variable (called the location variable).

}

\source{
See website from the refered article. The original data with 7,128 transcripts can be downloaded as an .RData file from \url{http://www.math.univ-toulouse.fr/~lecao/package.html}
}

\references{
Pomeroy, S.L. and Tamayo, P. and Gaasenbeek, M. and Sturla, L.M. and Angelo, M. and McLaughlin, M.E. and Kim, J.Y.H. and Goumnerova, L.C. and Black, P.M. and Lau, C. and others (2002). Prediction of central nervous system embryonal tumour outcome based on gene expression. \emph{Nature}, \bold{415}, 6870, 436--442. 


Hunt, L. and Jorgensen, M. (1999). Mixture model clustering using the MULTIMIX program. \emph{Australian \& New Zealand Journal of Statistics}, \bold{41}, 2, 154--171. 
}

\keyword{datasets}

\name{kmeans.init}
\encoding{latin1}
\alias{kmeans.init}
\title{
Initialization of the parameters with K-means
}
\description{
Parameters in integrativeME are first initialized via K-means clustering algorithm.
}


\usage{kmeans.init(data.cont)}

\arguments{
  \item{data.cont}{gene expression data, where a small number of genes has been selected beforehand. See \code{\link{integrativeME}}.}
}

\details{The number of clusters in K-means is set by default to 2. K-means is used to initialize the parameters before applying the EM algorithm in integrativeME.}

\value{
  \item{prop.kmeans}{proportions.}
  \item{means.kmeans}{estimated means.}
  \item{var.kmeans}{estimated variance-covariance matrix.}
  \item{continue}{boolean variable, set to FALSE if the number of observations is not sufficient to determine 2 clusters.}
}

\references{
Hartigan, J. A. and Wong, M. A. (1979). A K-means clustering algorithm. \emph{Applied Statistics}, \bold{28}, 100-108.



}

\author{Kim-Anh Le Cao}

\seealso{\code{\link{integrativeME}}, \code{\link{MEfunctions}}, \code{\link[stats]{kmeans}}.}

%\examples{}



\name{MEfunctions}
\encoding{latin1}
\alias{MEindep}
\alias{MElogreg}
\alias{MEloc}
\alias{MEcont}
\alias{MEfunctions}

\title{Internal gating functions for integrative Mixture of Experts methodology}

\description{
Different internal gating functions (or models) are proposed within mixture of experts to integrate gene expression data and clinical data in a binary classification framework.
}


\usage{

MEcont(jcross, train, test, n, nv, ng, indclass, data.cont, prop.kmeans, means.kmeans, var.kmeans)

MEindep(jcross, train, test, n, nq, nv, ng, indclass, data.cat, data.cont, prop.kmeans, means.kmeans, var.kmeans)

MElogreg(jcross, train, test, n, nq, nv, ng, indclass, data.cat, data.cont, prop.kmeans)

MEloc(jcross, train, test, n, nq, nv, ng, indclass, data.cat, data.cont, prop.kmeans, means.kmeans, var.kmeans, loc.ind)
}

\arguments{
  \item{jcross}{which cross validation sample}
  \item{train}{training samples}
  \item{test}{test samples}
  \item{n}{number of observations or samples}
  \item{nq}{number of clinical variables}
  \item{nv}{number of genes, the genes should be selected beforehand, see \code{\link{integrativeME}}.}
  \item{ng}{number of experts, should be set to 2 for a binary classification problem.}
  \item{indclass}{number of samples of class 0.}
  \item{data.cat}{clinical data (categorical).}
  \item{data.cont}{gene expression data.}
  \item{prop.kmeans}{proportions, initialized with k-means, see also \code{\link{kmeans.init}}.}
  \item{means.kmeans}{means, initialized with k-means, see also \code{\link{kmeans.init}}.}
  \item{var.kmeans}{variance-covariance matrix, initialized with k-means, see also \code{\link{kmeans.init}}.}
  \item{loc.ind}{index of the location variable in the case of the \code{MEloc} model.}


}

\details{Given a training set and a test set, the parameters in integrativeME are learnt via the EM algorithm and then tested. All three geting functions are included in the main program \code{integrativeME}}

\value{
	\item{prop}{estimated proportions.}
	\item{w}{weighted variable vector in the expert networks function.}
	\item{loglik}{loglikelihood of the model after several iterations.}
	\item{mat.gum}{main output that is used in \code{integrativeME} to predict the class label of each tested observation.}

}

\references{
Le Cao et al. (2009), submitted.

Ng, S.K. and McLachlan, G.J. (2008). Expert Networks with Mixed Continuous and Categorical Feature Variables: a Location Modeling Approach. \emph{Machine Learning Research Progress}, ed. Hanna Peters and Mia Vogel, 1--14

Hunt, L. and Jorgensen, M. (1999). Mixture model clustering using the MULTIMIX program. \emph{Australian \& New Zealand Journal of Statistics}, \bold{41}, 2, 154--171. 

}

\author{Kim-Anh Le Cao}

\seealso{\code{\link{integrativeME}}, \code{\link{kmeans.init}}}

%\examples{}



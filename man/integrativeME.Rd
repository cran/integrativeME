\name{integrativeME}
\encoding{latin1}
\alias{integrativeME}
\title{
Integrative mixture of experts (ME) to combine clinical variables and microarray data in a binary classification framework.
}
\description{
We have implemented and further developed mixture of experts (ME) to combine clinical (categorical) data and microarray data. In order to create a hybrid signature, gene selection is performed in a first step with various approaches. The selected genes and the clinical data are then combined with mixture of experts methodology in a classification framework and the evaluation of the performance is performed via k-fold cross-validation.}


\usage{
integrativeME(
	data.cat,
	data.cont,
	type,
	select = c('RF', 'student', 'sPLS'),
	method = c('logreg', 'indep', 'loc', 'cont'),
	loc.ind =  NULL,
	keepX = 5,
	ng = 2,
	mode.sPLS = NULL,
	fold = 10,
	kcv = 1                    # number of 10-fold cv    
)
}

\arguments{
  \item{data.cat}{clinical categorical data. The number of clinical factors is usually very small. Data should be sorted according to the \code{type} class vector.}
  \item{data.cont}{microarray data. Data should be sorted according to the \code{type} class vector.}
  \item{type}{vector indicating the class of each observation or sample. The vector should be coded 0 and 1 and sorted.}
  \item{select}{gene selection method to be applied in the first step. RF = random forests (wrapper method), student = t-test (filter method) and sPLS = sparse PLS, see also \code{\link[mixOmics]{spls}} to select genes according to the clinical factors.}
  \item{method}{variant of ME to use to combine both types of variables, see also \code{\link{MEfunctions}}.}
  \item{loc.ind}{if \code{method = 'loc'}, then the index of the location variable is needed.}
  \item{keepX}{number of genes to select. Should be set to a small value.}
  \item{ng}{number of experts to use in ME, set by default to 2 (i.e. the number of classes)}
  \item{mode.sPLS}{sPLS mode, \code{'regression'} to select genes according to \code{type}, \code{'canonical'} to select genes according to the clinical variables.}
  \item{fold}{number of folds in the cross-validation, by default set to 10}
  \item{kcv}{numbers of runs of the cross-validation.}
}

\details{
Clinical variables should all be categorical and should not contain any 0 values.

The samples in \code{data.cat} and \code{data.cont} should be ordered according to class vector \code{type}.

method = 'cont' is for continuous data only, whereas 'logreg', 'indep' and 'loc' combine both continuous and categorical variables together.

If \code{method = 'logreg'}, some warning messages may appear:  'fitted probabilities numerically 0 or 1 occurred' for binomial GLMs, see Venables \& Ripley (2002, pp. 197-8).


}

\value{
	\item{mean.error}{classification error rate for each cross-validation run.}
	\item{mat.predcited}{predicted class of each observation and each cross-validation run.}
}

\references{
Le Cao et al. (2009), submitted.

Ng, S.K. and McLachlan, G.J. (2008). Expert Networks with Mixed Continuous and Categorical Feature Variables: a Location Modeling Approach. \emph{Machine Learning Research Progress}, ed. Hanna Peters and Mia Vogel, 1--14

Hunt, L. and Jorgensen, M. (1999). Mixture model clustering using the MULTIMIX program. \emph{Australian \& New Zealand Journal of Statistics}, \bold{41}, 2, 154--171. 

}

\author{Kim-Anh Le Cao}

\seealso{\code{\link{MEfunctions}}, \code{\link{kmeans.init}}}

\examples{
\dontrun{
data(prostate)
data.cont = prostate$cont
data.cat = prostate$indep
type = prostate$type
type #check that type is sorted
# gene selection with t-test and ME model with a logistic regression in the gatin network:
res = integrativeME(
	data.cat = data.cat,
	data.cont = data.cont,
	type = type,
	select = 'student',
	method = 'logreg',
	keepX = 5,
	ng = 2,
	fold = 10,
	kcv = 1                    
)
}


}



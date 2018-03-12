\name{CalculatePAContinuous}
\alias{CalculatePAContinuous}

\title{Parallel Analysis for continuous data.}

\description{Obtains a parallel analysis for continuous data.}

\usage{
CalculatePAContinuous(dataMatrix, percentiles = 0.99, nReplicates = 200,
                      use = "complete.obs", algorithm = "pearson")
}


\arguments{
  \item{dataMatrix}{\code{matrix} or \code{data.frame} of binary or dichotomous variables.}
  \item{percentiles}{vector of percentiles to report.}
  \item{nReplicates}{number of simulations to produce for estimating the eigenvalues distribution under independence.}
  \item{use}{Missing value handling method: If \code{"complete.obs"}, remove observations with any missing data; if \code{"pairwise.complete.obs"}, compute each correlation using all observations with valid data for that pair of variables.}
  \item{algorithm}{string specifying the estimation algorithm. In the case of continuous variables, only the Pearson
    correlations are used. Ignored if different to \code{"pearson"}.}
}
\value{
  Returns a \code{list} object with the following:
  \item{observed}{\code{data.frame} containing the observed eigenvalues.}
  \item{percentiles}{\code{data.frame} containing the estimated percentiles of the eigenvalues distribution under
    independence.}
  \item{simulatedEigenValues}{\code{data.frame} containing the simulated eigenvalues under independence.}
 
}

\author{Carlos A. Arias \email{carias@icfes.gov.co} and Victor H. Cervantes \email{vcervantes@icfes.gov.co}}

\note{This is an auxiliary function for the \code{"PA"} function.}

\seealso{\code{\link{CalculatePABinary}}, \code{\link{CalculatePAOrdered}}, \code{\link{CalculatePAMixed}}, \code{\link{PA}}}

\examples{

# # Run Parallel analyis of numeric data (Iris)
data(iris)
continuousPA <- PA(iris[, -5], percentiles = c(0.90, 0.99), nReplicates = 200,
                   type = "continuous", algorithm = "pearson")
print(continuousPA)
}

\keyword{continuous}
\keyword{PA}

\encoding{latin1}

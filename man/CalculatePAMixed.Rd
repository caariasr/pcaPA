\name{CalculatePAMixed}
\alias{CalculatePAMixed}

\title{Parallel Analysis for numeric and ordered mixed data.}

\description{Obtains a parallel analysis for numeric and ordered mixed data.}

\usage{
CalculatePAMixed(dataMatrix, percentiles = 0.99, nReplicates = 200,
                 use = "complete.obs", algorithm = "polycor")
}

\arguments{
  \item{dataMatrix}{\code{matrix} or \code{data.frame} of ordered variables.}
  \item{percentiles}{vector of percentiles to report.}
  \item{nReplicates}{number of simulations to produce for estimating the eigenvalues distribution under independence.}
  \item{use}{Missing values handling method: If \code{"complete.obs"}, remove observations with any missing data; if \code{"pairwise.complete.obs"}, compute each correlation using all observations with valid data for that pair of variables.}
  \item{algorithm}{string specifying the correlation estimation method. In the case of mixed variables, only the polycor
    package is currently used, so this value must always be "polycor".}
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

\seealso{\code{\link{CalculatePABinary}}, \code{\link{CalculatePAContinuous}}, \code{\link{CalculatePAOrdered}}, \code{\link{PA}}}

\examples{
# # NOT RUN
# # Run Parallel analysis of mixed ordered and continuous data
# data(mixedScience)
# mixedPA <- PA(mixedScience, percentiles = c(0.90, 0.99),
#               nReplicates = 200, type = "mixed")
# print(mixedPA)
}

\keyword{mixed}
\keyword{PA}

\encoding{latin1}

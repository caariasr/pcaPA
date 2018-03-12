\name{CalculatePABinary}
\alias{CalculatePABinary}

\title{Parallel Analysis for Dichotomous Data.}

\description{Obtains a parallel analysis for dichotomous data.}

\usage{
CalculatePABinary(dataMatrix, percentiles = 0.99,nReplicates = 200,
                  use = "complete.obs", algorithm = "polycor")
}


\arguments{
  \item{dataMatrix}{\code{matrix} or \code{data.frame} of binary or dichotomous variables.}
  \item{percentiles}{vector of percentiles to report.}
  \item{nReplicates}{number of simulations to produce for estimating the eigenvalues distribution under independence.}
  \item{use}{Missing value handling method: If \code{"complete.obs"}, remove observations with any missing data; if \code{"pairwise.complete.obs"}, compute each correlation using all observations with valid data for that pair of variables.}
  \item{algorithm}{string specifying the correlation estimation algorithm. Polychoric correlation estimation method: \code{"polycor"} for estimates using the polycor package,
    \code{"polychoric"} for estimates using the \code{C++} function \code{Cpolychoric}.}
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

\seealso{\code{\link{CalculatePAOrdered}}, \code{\link{CalculatePAContinuous}}, \code{\link{CalculatePAMixed}}, \code{\link{PA}}, \code{\link{quantile.PA}}}

\examples{
# # NOT RUN
# # Run Parallel Analysis for binary data conforming to the Rasch model
# # using the polycor package
# data(simRaschData)
# binaryRaschPA <- PA(simRaschData, percentiles = c(0.95, 0.99),
#                     nReplicates = 200, type = "binary")
# print(binaryRaschPA)

# # Run Parallel Analysis for binary data conforming to the Rasch model
# # using the Cpolychoric C++ function
data(simRaschData)
binaryRaschPA <- PA(simRaschData, percentiles = c(0.95, 0.99), nReplicates = 200,
                    type = "binary", algorithm = "polychoric")
print(binaryRaschPA)

# # NOT RUN
# # Run Parallel Analysis for binary data conforming to the 2PL model
# # using the polycor package
# data(sim2plData)
# binary2plPA <- PA(sim2plData, percentiles = c(0.95, 0.99), nReplicates = 200,
#                   type = "binary")
# print(binary2plPA)

# # Run Parallel Analysis for binary data conforming to the 2PL model
# # using the polychoric C++ function
data(sim2plData)
binary2plPA <- PA(sim2plData, percentiles = c(0.95, 0.99), nReplicates = 200,
                  type = "binary", algorithm = "polychoric")
print(binary2plPA)
}

\keyword{dichotomous}
\keyword{PA}


\encoding{latin1}

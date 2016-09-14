\name{CalculatePAOrdered}
\alias{CalculatePAOrdered}

\title{Parallel Analysis for Ordered Data.}

\description{Obtains a parallel analysis for ordered data.}

\usage{
CalculatePAOrdered(dataMatrix, percentiles = 0.99, nReplicates = 200,
                   use = "complete.obs", algorithm = "polycor")
}

\arguments{
  \item{dataMatrix}{\code{matrix} or \code{data.frame} of ordered variables.}
  \item{percentiles}{vector of percentiles to report.}
  \item{nReplicates}{number of simulations to produce for estimating the eigenvalues distribution under independence.}
  \item{use}{Missing value handling method: If \code{"complete.obs"}, remove observations with any missing data; if \code{"pairwise.complete.obs"}, compute each correlation using all observations with valid data for that pair of variables.}
 \item{algorithm}{string specifying the correlation estimation algorithm. Polychoric correlation estimation method: \code{"polycor"} for estimates using the polycor package, \code{"polychoric"} for estimates using the \code{C++} function \code{Cpolychoric}.}
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

\seealso{\code{\link{CalculatePABinary}}, \code{\link{CalculatePAContinuous}}, \code{\link{CalculatePAMixed}}, \code{\link{PA}}}

\examples{
# # NOT RUN
# # Run Parallel analysis for ordered polytomous data using the polycor package
# data(Science)
# Science[, ] <- lapply(Science, as.ordered)
# orderedPA <- PA(Science, percentiles = c(0.90, 0.99), nReplicates = 200,
#                 type = "ordered")
# print(orderedPA)

# # Run Parallel analysis for ordered polytomous data using 
# # he polychoric C++ function
data(Science)
Science[, ] <- lapply(Science, as.ordered)
orderedPA   <- PA(Science, percentiles = c(0.90, 0.99), nReplicates = 200,
                 type = "ordered", algorithm = "polychoric")
print(orderedPA)


}

\keyword{ordered}
\keyword{PA}

\encoding{latin1}

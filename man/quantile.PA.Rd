\name{quantile.PA}
\alias{quantile.PA}

\title{Generate new quantiles based on given percentiles for a PA object.}

\description{Calculates arbitrary quantiles from the \code{simulatedEigenValues} in a \code{"PA"} object.}

\usage{\method{quantile}{PA}(x, percentiles = .99, ...)}

\arguments{
  \item{x}{An object of class \code{"PA"}.}
  \item{percentiles}{The new percentiles to replace the ones found in the \code{"PA"} object.}
  \item{\ldots}{Not used.}
}

\value{
  An object of class \code{"PA"} with the percentiles \code{data.frame} replaced by those from the percentiles
  argument.
}


\author{Carlos A. Arias \email{carias@icfes.gov.co} and Victor H. Cervantes \email{vcervantes@icfes.gov.co}} 

\seealso{\code{\link{PA}}, \code{\link{print.PA}}, \code{\link{Check.PA}}, \code{\link{CountEigen.PA}}, \code{\link{coef.PA}}, \code{\link{plot.PA}}}


\examples{
# # Using the polychoric C++ function
data(Science)
Science[, ] <- lapply(Science, as.ordered)
orderedPA   <- PA(Science, percentiles = c(0.90, 0.99), nReplicates = 200,
                  type = "ordered", algorithm = "polychoric")

print(orderedPA)  # Shows the eigenvalues as a matrix
head(orderedPA$simulatedEigenValues)  # Shows the first six iterations

# # Get different quantiles from a PA object
quantile(orderedPA, percentiles = c(0.1, 0.2, 0.5))
}

\keyword{quantile}
\keyword{PA}

\encoding{latin1}

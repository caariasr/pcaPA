\name{coef.PA}
\alias{coef.PA}

\title{Eigenvalue and percentile extraction of a \code{"PA"} object.}

\description{
  \code{coef} method for objects of class \code{"PA"}, produced by \code{PA}.}

\usage{
\method{coef}{PA}(object, ...)
}

\arguments{
  \item{object}{an object of class \code{"PA"}}
  \item{\ldots}{not used}
}

\value{
  An object of class \code{"matrix"} with the observed eigenvalues and the percentiles.
}

\author{Carlos A. Arias \email{carias@icfes.gov.co} and Victor H. Cervantes \email{vcervantes@icfes.gov.co}}

\seealso{\code{\link{PA}}, \code{\link{print.PA}}, \code{\link{Check.PA}}, \code{\link{CountEigen.PA}}, \code{\link{plot.PA}}, \code{\link{quantile.PA}}}

\examples{
# # Run Parallel Analysis for binary data conforming to the Rasch model
data(simRaschData)
binaryRaschPA <- PA(simRaschData, percentiles = c(0.95, 0.99), nReplicates = 200,
                    type = "binary", algorithm = "polychoric")
print(binaryRaschPA)
binaryRaschPAEigenValues <- coef(binaryRaschPA)  # Save the matrix of observed
                                                 # eigenvalues and estimated
                                                 # eigenvalue percentiles
binaryRaschPAEigenValues

}

\keyword{coef}
\keyword{methods}

\encoding{latin1}

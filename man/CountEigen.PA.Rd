\name{CountEigen.PA}
\alias{CountEigen.PA}

\title{Number of observed eigenvalues that exceed a given set of percentiles.}

\description{Counts the number of observed eigenvalues that exceed the given percentiles.}

\usage{CountEigen.PA(PA, percentiles = NULL)}

\arguments{
  \item{PA}{an object of class \code{"PA"}.}
  \item{percentiles}{the percentiles that ought to be plotted, defaults to those in the object.}
}

\value{
 A named numeric vector indicating the number of eigenvalues that are greater than the eigenvalues distribution
 percentiles under independence.
}


\author{Carlos A. Arias \email{carias@icfes.gov.co} and Victor H. Cervantes \email{vcervantes@icfes.gov.co}}

\seealso{\code{\link{PA}}, \code{\link{print.PA}}, \code{\link{coef.PA}}, \code{\link{Check.PA}}, \code{\link{plot.PA}}, \code{\link{quantile.PA}}}

\examples{

# # Run Parallel Analysis for binary data conforming to the Rasch model
data(simRaschData)
binaryRaschPA <- PA(simRaschData, percentiles = c(0.95, 0.99), nReplicates = 200,
                    type = "binary", algorithm = "polychoric")
print(binaryRaschPA)

# # Number of retained factors
nComponents <- CountEigen.PA(binaryRaschPA, percentiles = .99)
nComponents["p99"]

# # Run Parallel Analysis for binary data conforming to the 2PL model
data(sim2plData)
binary2plPA <- PA(sim2plData, percentiles = c(0.95, 0.99), nReplicates = 200,
                  type = "binary", algorithm = "polychoric")
print(binary2plPA)

# # Number of retained factors
nComponents <- CountEigen.PA(binary2plPA, percentiles = .99)
nComponents["p99"]

}

\keyword{PA}

\encoding{latin1}

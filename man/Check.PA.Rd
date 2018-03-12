\name{Check.PA}
\alias{Check.PA}

\title{Verifies that an object belongs to the \code{"PA"} class.}

\description{Checks if an object is of class \code{"PA"} (Parallel analysis).}

\usage{Check.PA(PA)}

\arguments{
  \item{PA}{An object to be checked for class \code{"PA"}.}
}

\value{
  Returns \code{TRUE} if input object if of class \code{"PA"} and has all the necessary components. Returns \code{FALSE} otherwise.
}

\author{Carlos A. Arias \email{carias@icfes.gov.co} and Victor H. Cervantes \email{vcervantes@icfes.gov.co}}

\seealso{\code{\link{PA}}, \code{\link{print.PA}}, \code{\link{coef.PA}}, \code{\link{CountEigen.PA}}, \code{\link{plot.PA}}, \code{\link{quantile.PA}}}


\examples{

# # Run Parallel Analysis for binary data conforming to the Rasch model
data(simRaschData)
binaryRaschPA <- PA(simRaschData, percentiles = c(0.95, 0.99), nReplicates = 200,
                    type = "binary", algorithm = "polychoric")
print(binaryRaschPA)

# # Check if binaryRaschPA is a PA object
Check.PA(binaryRaschPA)  # Should return TRUE
Check.PA(simRaschData)  # Should return FALSE

# # Run Parallel Analysis for binary data conforming to the 2PL model
data(sim2plData)
binary2plPA <- PA(sim2plData, percentiles = c(0.95, 0.99), nReplicates = 200,
                  type = "binary", algorithm = "polychoric")
print(binary2plPA)

# # Check if binary2plPA is a PA object
Check.PA(binary2plPA)  # Should return TRUE
Check.PA(simRaschData)  # Should return FALSE

}

\keyword{PA}

\encoding{latin1}

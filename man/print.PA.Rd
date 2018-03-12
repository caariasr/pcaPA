\name{print.PA}
\alias{print.PA}

\title{Print method for PA objects.}

\description{
  \code{print} method for objects of class \code{"PA"}. Prints the observed eigenvalues and the percentiles as a matrix.
}

\usage{\method{print}{PA}(x, digits = max(3, getOption("digits") - 3), observed = "Observed",
      percentile = "th percentile", position = "after", sep = "", ...)}

\arguments{
  \item{x}{an object of class \code{"PA"}.}
  \item{digits}{number of digits to print.}
  \item{observed}{Label for the observed data, default is "observed"}
  \item{percentile}{Graph title instead of default.}
  \item{position}{Position for the percentile label. \code{"after"} will position the label after the percentile number. \code{"before"} will position the label before the percentile number}
  \item{sep}{Character string to separate the label from the percentiles number.}
  \item{\ldots}{Not used.}
}

\value{
  NULL
}

\author{Carlos A. Arias \email{carias@icfes.gov.co} and Victor H. Cervantes \email{vcervantes@icfes.gov.co}}

\seealso{\code{\link{PA}}, \code{\link{plot.PA}}, \code{\link{Check.PA}}, \code{\link{CountEigen.PA}}, \code{\link{coef.PA}}, \code{\link{quantile.PA}}}

\examples{
# # Run Parallel Analysis for binary data conforming to the Rasch model
# # using the polycor package
data(simRaschData)
binaryRaschPA <- PA(simRaschData, percentiles = c(0.95, 0.99), nReplicates = 200,
                    type = "binary", algorithm = "polychoric")
print(binaryRaschPA)
}

\keyword{print}
\keyword{methods}

\encoding{latin1}

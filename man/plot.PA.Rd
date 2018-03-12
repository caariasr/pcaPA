\name{plot.PA}
\alias{plot.PA}

\title{Plot method for PA objects.}

\description{\code{plot} method for objects of class \code{PA}. Plots the scree plot for a \code{"PA"} object for the
  selected percentiles using \code{ggplot}.
}


\usage{
\method{plot}{PA}(x, percentiles = NULL, main = NULL, xlab = NULL, ylab = NULL,
      groupLabel = NULL, colour = TRUE, linetype = TRUE, observed = "Observed",
      percentile = "th percentile", position = "after", sep = "", ...)
}

\arguments{
  \item{x}{an object of class \code{"PA"}.}
  \item{percentiles}{The percentiles that ought to be plotted. Defaults to those in the PA object.}
  \item{main}{Graph title instead of default.}
  \item{xlab}{Label for \code{x} axis instead of default.}
  \item{ylab}{Label for \code{y} axis instead of default.}
  \item{groupLabel}{Legend box name instead of default.}
  \item{colour}{Logical indicating whether to identify the observed eigenvalues and percentiles by colour.}
  \item{linetype}{Logical indicating whether to identify the observed eigenvalues and percentiles by linetype.}
  \item{observed}{Label for the observed data, default is "observed"}
  \item{percentile}{Graph title instead of default.}
  \item{position}{Position for the percentile label. \code{"after"} will position the label after the percentile number. \code{"before"} will position the label before the percentile number}
  \item{sep}{Character string to separate the label from the percentiles number.}
  \item{\ldots}{Not used.}}

\value{
  \code{ggplot} object for plotting the scree plot.
}

\author{Carlos A. Arias \email{carias@icfes.gov.co} and Victor H. Cervantes \email{vcervantes@icfes.gov.co}} 

\seealso{\code{\link{PA}}, \code{\link{print.PA}}, \code{\link{Check.PA}}, \code{\link{CountEigen.PA}}, \code{\link{coef.PA}}, \code{\link{quantile.PA}}}


\examples{
# # Run Parallel Analysis for binary data conforming to the Rasch model
# # using the polychoric C++ function
data(simRaschData)
binaryRaschPA <- PA(simRaschData, percentiles = c(0.95, 0.99), nReplicates = 200,
                    type = "binary", algorithm = "polychoric")
print(binaryRaschPA)
plot(binaryRaschPA, percentiles = 0.99, groupLabel = "") # Plots the scree-plot
                                                         # with the 99th percentile
}

\keyword{plot}
\keyword{PA}

\encoding{latin1}

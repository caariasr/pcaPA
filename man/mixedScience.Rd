\name{mixedScience}
\alias{mixedScience}
\docType{data}
\title{
Simulated data from a normal distribution added to the Science data set from package "ltm".
}
\description{
Mixed ordered (Science dataset from package "ltm" coming from the Consumer Protection and Perceptions of Science and
Technology section of the 1992 Euro-Barometer Survey) and continuous (simulated from a normal distribution) data.
}
\usage{data(mixedScience)}
\format{
  A data frame with 392 observations on the following 12 variables:
  \describe{
    \item{\code{Comfort}}{an ordered factor with levels \code{strongly disagree} < \code{disagree} < \code{agree} < \code{strongly agree}}
    \item{\code{Environment}}{an ordered factor with levels \code{strongly disagree} < \code{disagree} < \code{agree} < \code{strongly agree}}
    \item{\code{Work}}{an ordered factor with levels \code{strongly disagree} < \code{disagree} < \code{agree} < \code{strongly agree}}
    \item{\code{Future}}{an ordered factor with levels \code{strongly disagree} < \code{disagree} < \code{agree} < \code{strongly agree}}
    \item{\code{Technology}}{an ordered factor with levels \code{strongly disagree} < \code{disagree} < \code{agree} < \code{strongly agree}}
    \item{\code{Industry}}{an ordered factor with levels \code{strongly disagree} < \code{disagree} < \code{agree} < \code{strongly agree}}
    \item{\code{Benefit}}{an ordered factor with levels \code{strongly disagree} < \code{disagree} < \code{agree} < \code{strongly agree}}
    \item{\code{X1}}{a numeric vector}
    \item{\code{X2}}{a numeric vector}
    \item{\code{X3}}{a numeric vector}
    \item{\code{X4}}{a numeric vector}
    \item{\code{X5}}{a numeric vector}
  }
}

\examples{
# # NOT RUN
# data(mixedScience)
# mixedPA <- PA(mixedScience, percentiles = c(0.90, 0.99), nReplicates = 200,
#               type = "mixed")
# print(mixedPA)
}
\keyword{datasets}

\references{
   Karlheinz, R. and Melich, A. (1992). Euro-Barometer 38.1:   Consumer Protection and Perceptions of Science and Technology. INRA (Europe), Brussels. [computer file]

Rizopoulos, D. (2006). ltm: An R package for Latent Variable Modelling and Item Response Theory Analyses. Journal of
Statistical Software, 17(5), 1--25.
}
\encoding{latin1}

\name{sim2plData}
\alias{sim2plData}
\docType{data}
\title{
Simulated data conforming to the 2pl model.
}
\description{
Dichotomous data simulated from a 2pl model.
}
\usage{data(sim2plData)}
\format{
  The format is:
 num [1:300, 1:15] 1 0 1 1 1 1 0 1 1 0 ...
}

\source{
Data simulated using the eRm package.
}

\examples{
data(sim2plData)
binary2plPA <- PA(sim2plData, percentiles = c(0.95, 0.99), nReplicates = 200,
                  type = "binary", algorithm = "polychoric")
print(binary2plPA)
}
\keyword{datasets}

\encoding{latin1}

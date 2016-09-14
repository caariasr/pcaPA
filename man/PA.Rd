\name{PA}
\alias{PA}


\title{General function to perform parallel analysis of continuous, ordered or mixed type data.}

\description{Creates an object of class PA using one of the PA functions. It performs a parallel analysis of continuous,
orderd (including dichotomous/binary as a special case) or mixed type of data associated with a principal components
analysis. Polychoric correlations among ordered variables, Pearson correlations among continuous variables and
polyserial correlation between mixed type variables (one ordered and one continuous) are used. Whenever the use of
polyserial or polychoric correlations yields a non positive definite correlation matrix, the resulting matrix is
transformed into the nearest positive definite matrix by \code{nearcor}.}

\usage{
PA(dataMatrix, percentiles = 0.99, nReplicates = 200, type = "continuous",
   use = "complete.obs", algorithm = "polycor")
}

\arguments{
  \item{dataMatrix}{\code{matrix} or \code{data.frame} of continuous numeric variables.}
  \item{percentiles}{vector of percentiles to report.}
  \item{nReplicates}{number of simulations to produce for estimating the eigenvalues distribution under independence.}
  \item{type}{Data type: \code{"continuous"} if the data is continuous. \code{"binary"} if the data is dichotomous. \code{"ordered"} if the data is ordinal. \code{"mixed"} if the data is mixed rodered and continuous.} 
   \item{use}{Missing value handling method: If \code{"complete.obs"}, remove observations with any missing data; if \code{"pairwise.complete.obs"}, compute each correlation using all observations with valid data for that pair of variables.}
  \item{algorithm}{string specifying the correlation estimation algorithm. Polychoric correlation estimation method:
    \code{"polycor"} for estimates using the polycor package, \code{"polychoric"} for estimates using the \code{C++}
    function \code{polychoric}. Pearson correlation estimation: \code{"pearson"}, only used for continuous data.}
}


\details{
This function generates an object of (S3) class \code{"PA"} by using one of the four support functions
(\code{\link{CalculatePAContinuous}}, \code{\link{CalculatePAOrdered}}, \code{\link{CalculatePABinary}},
\code{\link{CalculatePAMixed}}). As noted by Presaghi & Desimoni (2011), when using the parallel analysis approach to
select the number of components to retain, observed and simulated correlation matrices ought to be of the same kind.
The purporse with this package is to provide a more flexible framework to work with parallel anaylisis data.
For the case where all variables are ordered polytomous the \code{C++} function \code{polychoric} 
or the function \code{"polychor"} from package \code{"polycor"} may be used.
Pearson correlation from \code{"cor"} function is used; finally, for sets that mix variables of both types, the \code{"polycor"} package is 
used. Furthermore, objects of class \code{"PA"} retain the eigenvalues for simulated correlation matrices allowing the
user to avoid multiple simulations for changing certain parameters such as the specific qunatile used to decide on the
retained number of components; also some missing observations may be handled by the \code{"use"} parameter.
}


\value{
  An object of class \code{"PA"} with the following:
  \item{observed}{\code{data.frame} containing the observed eigenvalues.}
  \item{percentiles}{\code{data.frame} containing the estimated percentiles of the eigenvalues distribution under
    independence.}
  \item{simulatedEigenValues}{\code{data.frame} containing the simulated eigenvalues under independence.}
}


\author{Carlos A. Arias \email{carias@icfes.gov.co} and Victor H. Cervantes \email{vcervantes@icfes.gov.co}}

\note{The algorithm \code{"polychoric"} is only implemented for binary/dichotomous and ordered polytomous data. When
  mixed type data are used only \code{"polycor"} is currently available.}

\references{
Horn, J. L. (1965). A rationale and test for the number of factors in factor analysis. Psychometrika, 30, 179--185.

Glorfeld, L. W. (1995). An Improvement on Horn's Parallel Analysis Methodology for Selecting the Correct Number of Factors to Retain. Educational and Psychological Measurement, 55(3), 377--393.
}

\seealso{\code{\link{CalculatePAContinuous}}, \code{\link{CalculatePAOrdered}}, \code{\link{CalculatePABinary}},
  \code{\link{CalculatePAMixed}}, \code{\link{print.PA}}, \code{\link{plot.PA}}, \code{\link{coef.PA}},
  \code{\link{quantile.PA}}, \code{\link{CountEigen.PA}}, \code{\link{Check.PA}}}

\examples{
# # NOT RUN
# # Run Parallel Analysis for binary data conforming to the Rasch model
# # using the polycor package
# data(simRaschData)
# binaryRaschPA <- PA(simRaschData, percentiles = c(0.95, 0.99), nReplicates = 200,
#                     type = "binary")
# print(binaryRaschPA)

# # Run Parallel Analysis for binary data conforming to the Rasch model
# # using the polychoric C++ function
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

# # NOT RUN
# # Run Parallel analysis for ordered polytomous data using the polycor package
# data(Science)
# Science[, ] <- lapply(Science, as.ordered)
# orderedPA <- PA(Science, percentiles = c(0.90, 0.99), nReplicates = 200,
#                 type = "ordered")
# print(orderedPA)

# # Run Parallel analysis for ordered polytomous data using the polychoric C++ function
data(Science)
Science[, ] <- lapply(Science, as.ordered)
orderedPA   <- PA(Science, percentiles = c(0.90, 0.99), nReplicates = 200,
                  type = "ordered", algorithm = "polychoric")
print(orderedPA)

# # NOT RUN
# # Run Parallel analysis of mixed ordered and continuous data
# data(mixedScience)
# mixedPA <- PA(mixedScience, percentiles = c(0.90, 0.99), nReplicates = 200,
#               type = "mixed")
# print(mixedPA)

# # Run Parallel analyis of numeric data (Iris)
data(iris)
continuousPA <- PA(iris[, -5], percentiles = c(0.90, 0.99), nReplicates = 200,
                   type = "continuous", algorithm = "pearson")
print(continuousPA)
}

\keyword{Parallel Analysis}
\keyword{ParallelA}
\keyword{PA}

\encoding{latin1}

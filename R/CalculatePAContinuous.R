################################################################################
# # Function for parallel analysis of correlation matrices of Continuous data
################################################################################
CalculatePAContinuous <- function (dataMatrix, percentiles = 0.99, nReplicates = 200, use = "complete.obs",
                                   algorithm = "pearson") {
  # # Obtains a parallel analysis for continuous data
  # #
  # # Arg:
  # #  dataMatrix: matrix or data.frame of continuous numeric variables
  # #  percentiles: vector of percentiles to report
  # #  nReplicates: number of simulations to produce for estimating the eigenvalues distribution under independence
  # #  use: Missing value handling method: If "complete.obs", remove observations with any missing data; if "pairwise.complete.obs", compute each correlation using all observations with valid data for that pair of variables.
  # #  algorithm: string specifying the matrix for which the PCA analysis is conducted. Covariance matrix: "cov", or correlation matrix: "pearson".
  # #
  # # Ret:
  # #  parallelList: a list with data.frames observed, percentiles and simulated data.
  # #   observed: data.frame containing the observed eigenvalues
  # #   percentiles: data.frame containing the estimated percentiles of the eigenvalues distribution under independence
  # #   simulatedEigenValues: data.frame containing the simulated eigenvalues under independence


  ################################################################################
  # # Data verification
  ################################################################################
  if (!is.matrix(dataMatrix) & !is.data.frame(dataMatrix)) {
    stop("dataMatrix must be a matrix or a data.frame")
  }
  if (!is.data.frame(dataMatrix)) {
    dataMatrix <- data.frame(dataMatrix)
  }

  isNumericData <- all(sapply(dataMatrix, is.numeric))
  if (!isNumericData) {
    stop("All variables in dataMatrix must be numeric")
  }

  if (use %in% c("everything", "all.obs")) {
    remUse <- FALSE
  } else {
    remUse <- TRUE
  }

  if ((algorithm != "cov") & (algorithm != "pearson")) {
    warning("Only covariances or Pearson correlations are used for continuous data.\nUsing Pearson correlation.")
  }

  ################################################################################
  # # Data information
  ################################################################################
  nObservations  <- nrow(dataMatrix)
  nVariables     <- ncol(dataMatrix)
  dataMeans      <- colMeans(dataMatrix, na.rm = remUse)
  dataSds        <- apply(dataMatrix, 2, sd, na.rm = remUse)
  if (algorithm == 'cov') {
    datCorrelation <- cov(dataMatrix, use = use)
  } else {
    datCorrelation <- cor(dataMatrix, use = use)
  }
  datEigenValues <- eigen(datCorrelation)$values

  observed <- data.frame(orderEigenValues = 1:nVariables,
                         typeEigenValues  = "Observed",
                         eigenValues      = datEigenValues,
                         stringsAsFactors = TRUE)

  ################################################################################
  # # Simulate correlation matrices under independence
  ################################################################################
  simulatedEigenValues           <- matrix(nrow = nReplicates, ncol = nVariables)
  rownames(simulatedEigenValues) <- paste("iter", 1:nReplicates, sep = "")
  colnames(simulatedEigenValues) <- 1:nVariables

  for (ii in 1:nReplicates) {
    simulatedData <- matrix(rnorm(nObservations * nVariables,
                                  mean = dataMeans, sd = dataSds),
                            ncol = nVariables, nrow = nObservations)

    simulatedEigenValues[ii, ] <- eigen(cor(simulatedData, use = use))$values
  }
  simulatedEigenValues <- data.frame(simulatedEigenValues)

  ################################################################################
  # # Obtain percentiles of simulated data
  ################################################################################
  estimatedPercentiles <- sapply(simulatedEigenValues, quantile, percentiles)

  if (length(percentiles) == 1) {
    estimatedPercentiles <- data.frame(orderEigenValues = 1:nVariables,
                                       typeEigenValues  = 100 * percentiles,
                                       eigenValues      = estimatedPercentiles)
    rownames(estimatedPercentiles) <- 1:nVariables
  } else {
    estimatedPercentiles <- data.frame(t(estimatedPercentiles))
    pValues <- grep("^X\\d+\\.$", names(estimatedPercentiles))
    names(estimatedPercentiles)[pValues] <- gsub("^X(\\d+)\\.$", "p.\\1", names(estimatedPercentiles)[pValues])
    estimatedPercentiles <- reshape(estimatedPercentiles, direction = "long",
                                    varying = seq(along = names(estimatedPercentiles)))
    estimatedPercentiles <- data.frame(orderEigenValues = estimatedPercentiles[, "id"],
                                       typeEigenValues  = estimatedPercentiles[, "time"],
                                       eigenValues      = estimatedPercentiles[, "p"])
  }

  ################################################################################
  # # Output
  ################################################################################
  parallelList <- list(observed = observed, percentiles = estimatedPercentiles,
                       simulatedEigenValues = simulatedEigenValues)

  return(parallelList)
}

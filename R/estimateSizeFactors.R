#' Prepare Count Data for Differential Expression Analysis
#'
#' Converts various supported input types to a standardized list format for downstream
#' differential expression analysis. Supports `matrix`, `data.frame`, `DGEList`,
#' `DESeqDataSet`, and `SummarizedExperiment` objects.
#'
#' @param data A numeric matrix, data.frame, or supported object containing counts.
#' @param group A vector of group labels for the columns/samples of `data`. Must be the same length as the number of columns in `data`.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{counts}{An integer matrix of counts.}
#'   \item{samples}{A data.frame containing sample-level metadata: `group`, `lib.size`, and `size.factor`.}
#' }
#'
#' @details
#' This function performs input validation.
#' - Checks for non-negative numeric values and absence of NA.
#' - Ensures group labels match the number of samples.
#' - Automatically assigns column names if missing.
#' - Returns a list suitable for use with hurdle model-based DE functions.
#'
#' @examples
#' \dontrun{
#' # Example with a matrix
#' set.seed(123)
#' mat <- matrix(rnbinom(30, size = 5, mu = 5), nrow = 5)
#' grp <- c("A", "A", "A", "B", "B", "B")
#' y <- prepareDGE(mat, grp)
#'
#' # Example with a SummarizedExperiment
#' # se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
#' # y <- prepareDGE(se, grp)
#' }
#'
#' @export
prepareDGE <- function(data, group) {
  # 1. convert supported input types to a count matrix
  if (inherits(data, "DGEList")) {
    data <- data$counts

  } else if (inherits(data, "DESeqDataSet") || inherits(data, "SummarizedExperiment")) {
    # Check if SummarizedExperiment is installed
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop(
        "Package 'SummarizedExperiment' is required to process this input type. ",
        "Please install it or provide a matrix/data.frame instead."
      )
    }
    data <- SummarizedExperiment::assay(data)

  } else if (is.data.frame(data)) {
    data <- as.matrix(data)

  } else if (!is.matrix(data)) {
    stop(paste(
      "'data' must be one of: matrix, data.frame, DGEList,",
      "DESeqDataSet, or SummarizedExperiment."
    ))
  }

  # 2. basic checks
  if (is.null(colnames(data))) {
    colnames(data) <- paste0("sample", seq_len(ncol(data)))
  }

  if (!is.numeric(data)) {
    stop("'data' must contain numeric counts.")
  }

  if (length(group) != ncol(data)) {
    stop("'group' must have length equal to ncol(data).")
  }

  if (any(is.na(data))) {
    stop("'data' contains NA values.")
  }

  if (any(data < 0)) {
    stop("'data' must contain non-negative counts.")
  }

  # optional: enforce integer counts for DESeq2 compatibility
  if (any(abs(data - round(data)) > 1e-8)) {
    stop("'data' must contain integer counts.")
  }
  storage.mode(data) <- "integer"

  # 3. build samples exactly as before
  lib.size <- colSums(data)
  size.factor <- rep(1, ncol(data))
  samples <- data.frame(group = as.factor(group),
                        lib.size = lib.size,
                        size.factor = size.factor)
  rownames(samples) <- colnames(data)

  return(list(counts = data, samples = samples))
}


#' Estimate Size Factors for Normalization
#'
#' Computes sample-specific size factors for long-read RNA-Seq data, used to
#' normalize counts for differential expression analysis.
#'
#' @param y A count matrix (`matrix` or `data.frame`) or the output of \code{\link{prepareDGE}()}.
#'   If a matrix/data.frame is provided, the function assumes two equal-sized groups.
#' @param type Character string specifying the method for estimating size factors:
#'   \describe{
#'     \item{"poscounts"}{Geometric mean-based method.}
#'     \item{"ratio"}{Simple ratio method using the mean of log counts.}
#'   }
#'   Default is \code{"poscounts"}.
#' @param locfunc Function to summarize log-ratios across genes. Defaults to \code{\link[stats]{median}}.
#'
#' @return A list (same structure as \code{\link{prepareDGE}()} output) with:
#'   \describe{
#'     \item{counts}{Original count matrix (integer).}
#'     \item{samples}{Data frame with sample information, updated \code{size.factor}.}
#'     \item{baseMean}{Normalized mean of counts per gene.}
#'   }
#'
#' @details
#' This function implements two methods for size factor estimation:
#'
#' * **poscounts**: Computes a geometric mean of positive counts per gene,
#'   then calculates ratios for each sample. Normalizes so that the geometric mean of size factors equals 1.
#' * **ratio**: Uses the mean of log-counts per gene across samples to compute ratios.
#'
#' The function automatically normalizes counts using the estimated size factors
#' and stores gene-level normalized means in \code{baseMean}.
#'
#' @examples
#' \dontrun{
#' # Using a count matrix
#' #' set.seed(123)
#' mat <- matrix(rnbinom(30, size = 5, mu = 5), nrow = 5)
#' grp <- c("A", "A", "A", "B", "B", "B")
#' y <- prepareDGE(mat, grp)
#' y <- sizeFactorsEst(y, type = "poscounts")
#' }
#'
#' @importFrom stats median
#' @export
sizeFactorsEst <- function(y, type = c("poscounts", "ratio"),
                           locfunc = stats::median) {
  type <- match.arg(type)

  if (is.data.frame(y) || is.matrix(y)) {
    if (ncol(y) %% 2 != 0) {
      stop("When 'y' is a matrix/data.frame, ncol(y) must be even to assume two equal-size groups.")
    }
    reps <- ncol(y) / 2
    groups <- factor(c(rep("Group1", reps), rep("Group2", reps)))
    message("Input is a matrix/data.frame, assuming two equal-size groups.")
    y <- prepareDGE(y, groups)
  }

  if (!is.list(y) || is.null(y$counts) || is.null(y$samples)) {
    stop("'y' must be either a matrix/data.frame or the output of prepareDGE().")
  }

  obs <- y$counts

  if (!is.matrix(obs)) {
    stop("y$counts must be a matrix.")
  }

  if (any(is.na(obs))) {
    stop("y$counts contains NA values.")
  }

  if (any(obs < 0)) {
    stop("y$counts must contain non-negative counts.")
  }

  if (any(colSums(obs) == 0)) {
    stop("at least one sample has total count 0, cannot estimate size factors.")
  }

  if (type == "poscounts") {
    geoMeanNZ <- function(x) {
      if (all(x == 0)) {
        0
      } else {
        exp(sum(log(x[x > 0])) / length(x))
      }
    }
    geoMeans <- apply(obs, 1, geoMeanNZ)
    loggeomeans <- log(geoMeans)

  } else if (type == "ratio") {
    loggeomeans <- rowMeans(log(obs))
  }

  if (all(is.infinite(loggeomeans))) {
    stop("every gene contains at least one zero, cannot compute log geometric means")
  }

  sf <- apply(obs, 2, function(cnts) {
    idx <- is.finite(loggeomeans) & cnts > 0
    if (!any(idx)) {
      NA_real_
    } else {
      exp(locfunc((log(cnts) - loggeomeans)[idx]))
    }
  })

  if (any(is.na(sf)) || any(sf <= 0)) {
    stop("failed to estimate valid size factors.")
  }

  if (type == "poscounts") {
    sf <- sf / exp(mean(log(sf)))
  }

  y$samples$size.factor <- sf

  norm_counts <- sweep(obs, 2, y$samples$size.factor, "/")
  y$baseMean <- rowMeans(norm_counts)

  return(y)
}

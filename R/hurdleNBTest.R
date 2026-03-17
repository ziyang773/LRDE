#' @importFrom stats optim pnorm pchisq plogis dnbinom dnorm qlogis
NULL
#' Negative Log-Likelihood for Hurdle NB with Fixed Zero Probabilities and Dispersion
#'
#' Internal function to compute the negative log-likelihood of a hurdle
#' negative binomial model where both zero probabilities (\code{theta}) and
#' dispersion (\code{phi}) are fixed. Only mean parameters are optimized.
#'
#' @param params Numeric vector:
#'   \describe{
#'     \item{1}{If \code{one_group = TRUE}, log(mu) shared across groups.}
#'     \item{1}{If \code{one_group = FALSE}, log(mu_A) for group A.}
#'     \item{2}{If \code{one_group = FALSE}, delta_rate = log(mu_B) - log(mu_A).}
#'   }
#' @param counts Numeric vector of observed counts for a single gene.
#' @param groups Factor indicating group membership for each observation.
#' @param offset Optional numeric vector of normalization factors. If \code{NULL},
#'   defaults to 1 for all observations.
#' @param fixed_theta_A Numeric, fixed non-zero probability for group A.
#' @param fixed_theta_B Numeric, fixed non-zero probability for group B.
#' @param fixed_phi Numeric, fixed NB dispersion parameter (size).
#' @param one_group Logical; if \code{TRUE}, a single mean is shared across both groups.
#'
#' @return Numeric, the negative log-likelihood value.
#'
#' @details
#' The likelihood consists of:
#' \itemize{
#'   \item Bernoulli component for zero vs non-zero counts using fixed theta
#'   \item Truncated negative binomial component for positive counts using fixed phi
#' }
#'
#' This function is used internally during optimization when only mean parameters
#' are estimated, while theta and phi are kept constant.
#'
#' @keywords internal
nll_hurdle_fixed_theta_phi <- function(params, counts, groups, offset = NULL,
                                       fixed_theta_A, fixed_theta_B,
                                       fixed_phi, one_group = TRUE) {
  eps <- 1e-12

  # --- 0. Setup ---
  if (is.null(offset)) offset <- rep(1, length(counts))
  levs <- levels(groups)

  # --- 1. Map Parameters ---
  # Only log_rates are in params now
  if (one_group) {
    # params = [log_rate_shared]
    log_rate_A <- params[1]
    log_rate_B <- params[1]
  } else {
    # params = [log_rate_A, delta_rate]
    log_rate_A <- params[1]
    delta_rate <- params[2]
    log_rate_B <- log_rate_A + delta_rate
  }

  # Map fixed thetas and rates to observations
  is_grp_A <- (groups == levs[1])
  thetas <- ifelse(is_grp_A, fixed_theta_A, fixed_theta_B)
  thetas <- pmin(pmax(thetas, eps), 1 - eps)
  rates  <- ifelse(is_grp_A, exp(log_rate_A), exp(log_rate_B))
  mus    <- rates * offset

  # --- 2. Likelihood Calculation ---
  is_zero <- counts == 0

  # A. Binary Part (using fixed thetas)
  ll_binary <- sum(log1p(-thetas[is_zero])) + sum(log(thetas[!is_zero]))

  # B. Truncated Negative Binomial Part
  pos_counts <- counts[!is_zero]
  pos_mus    <- mus[!is_zero]

  # Standard NB log-likelihood
  ll_nb <- sum(dnbinom(pos_counts, mu = pos_mus, size = fixed_phi, log = TRUE))

  # Truncation Adjustment: log(1 / (1 - P(Y=0)))
  prob_zero_nb     <- dnbinom(0, mu = pos_mus, size = fixed_phi)
  prob_zero_nb     <- pmin(pmax(prob_zero_nb, 0), 1 - eps)
  trunc_adjustment <- -sum(log1p(-prob_zero_nb))

  return(-(ll_binary + ll_nb + trunc_adjustment))
}


#' Likelihood Ratio Test for Hurdle Negative Binomial Model
#'
#' Performs gene-wise likelihood ratio tests (LRT) for differential expression
#' using a hurdle negative binomial model with fixed zero probabilities and
#' tag-wise dispersions.
#'
#' @param y A list-like object returned from \code{tagwiseEst()} containing:
#'   \describe{
#'     \item{counts}{Numeric matrix of gene expression counts (genes x samples).}
#'     \item{samples}{Data frame with columns \code{group} (factor) and \code{size.factor} (numeric).}
#'     \item{tagwise.disp}{Numeric vector of estimated tag-wise dispersions.}
#'     \item{zero_prob_matrix}{Matrix of zero probabilities per group per gene.}
#'   }
#'
#' @return The input object \code{y} with two additional elements:
#'   \describe{
#'     \item{lrt_stats}{Numeric vector of LRT statistics for each gene.}
#'     \item{p.values}{Numeric vector of corresponding p-values.}
#'   }
#'
#' @details
#' For each gene:
#' \itemize{
#'   \item The null model assumes a single shared mean across groups.
#'   \item The alternative model estimates group-specific means.
#'   \item Zero probabilities and dispersions are fixed from prior estimates.
#'   \item When one group has all zero counts, a one-sided Z-test is applied instead.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' mat <- matrix(rnbinom(30, size = 5, mu = 5), nrow = 5)
#' grp <- c("A", "A", "A", "B", "B", "B")
#' y <- prepareDGE(mat, grp)
#' y <- sizeFactorsEst(y)
#' y <- tagwiseEst(y)
#' y <- hurdle_LRT(y)
#' }
#'
#' @export
hurdle_LRT <- function(y){
  n_row <- nrow(y$counts)
  lrt_stats <- numeric(n_row)
  p_values  <- numeric(n_row)

  groups_factor <- factor(y$samples$group)
  offsets       <- y$samples$size.factor
  levs          <- levels(groups_factor)

  for (i in 1:n_row) {
    current_counts <- as.numeric(y$counts[i, ])
    current_phi <- 1 / y$tagwise.disp[i]

    # Retrieve the fixed priors for this gene
    f_theta_A <- 1 - y$zero_prob_matrix[i, "prob_A"]
    f_theta_B <- 1 - y$zero_prob_matrix[i, "prob_B"]

    all_zero_A <- all(current_counts[groups_factor == levs[1]] == 0)
    all_zero_B <- all(current_counts[groups_factor == levs[2]] == 0)

    # Initial estimate for log_rate (shared)
    obs_mean <- mean(current_counts[current_counts > 0] / offsets[current_counts > 0], na.rm=TRUE)
    if(is.na(obs_mean) || obs_mean <= 0) obs_mean <- 0.1

    tryCatch({
      # --- Model 1: Null (1 Shared Mean) ---
      # We estimate 1 parameter: log_rate
      null_model <- optim(
        par = log(obs_mean),
        fn = nll_hurdle_fixed_theta_phi,
        counts = current_counts,
        groups = groups_factor,
        offset = offsets,
        fixed_theta_A = f_theta_A,
        fixed_theta_B = f_theta_B,
        fixed_phi = current_phi,
        one_group = TRUE,
        method = "BFGS",
        hessian = TRUE
      )

      if (all_zero_A || all_zero_B) {
        # when one group has all zeros, test the mean of the other group
        se_beta <- sqrt(1 / null_model$hessian[1, 1])

        # Z-test: (Observed log-rate - "Zero" Floor) / SE
        # testing if the shared mean is significantly above a floor
        z_stat <- (null_model$par - 1) / se_beta

        lrt_stats[i] <- z_stat
        # One-sided test (upper tail)
        p_values[i]  <- pnorm(z_stat, lower.tail = FALSE)

      } else {
        # --- Model 2: Alt (Standard LRT) ---
        m_A <- mean(current_counts[groups_factor == levs[1] & current_counts > 0] /
                      offsets[groups_factor == levs[1] & current_counts > 0], na.rm=TRUE)
        m_B <- mean(current_counts[groups_factor == levs[2] & current_counts > 0] /
                      offsets[groups_factor == levs[2] & current_counts > 0], na.rm=TRUE)

        alt_model <- optim(
          par = c(log(m_A), log(m_B) - log(m_A)),
          fn = nll_hurdle_fixed_theta_phi,
          counts = current_counts,
          groups = groups_factor,
          offset = offsets,
          fixed_theta_A = f_theta_A,
          fixed_theta_B = f_theta_B,
          fixed_phi = current_phi,
          one_group = FALSE,
          method = "BFGS"
        )

        lrt <- 2 * (null_model$value - alt_model$value)
        lrt_stats[i] <- max(0, lrt)
        p_values[i]  <- pchisq(lrt_stats[i], df = 1, lower.tail = FALSE)
      }

    }, error = function(e) {
      lrt_stats[i] <<- NA
      p_values[i]  <<- NA
    })
  }

  y$lrt_stats <- lrt_stats
  y$p.values  <- p_values
  return(y)
}


#' Wald Test for Hurdle Negative Binomial Model
#'
#' Performs gene-wise Wald tests for differential expression using a hurdle
#' negative binomial model with fixed zero probabilities and tag-wise dispersions.
#'
#' @param y A list-like object returned from \code{tagwiseEst()} containing:
#'   \describe{
#'     \item{counts}{Numeric matrix of gene expression counts (genes x samples).}
#'     \item{samples}{Data frame with columns \code{group} (factor) and \code{size.factor} (numeric).}
#'     \item{tagwise.disp}{Numeric vector of estimated tag-wise dispersions.}
#'     \item{zero_prob_matrix}{Matrix of zero probabilities per group per gene.}
#'   }
#'
#' @return The input object \code{y} with two additional elements:
#'   \describe{
#'     \item{wald_stats}{Numeric vector of Wald statistics for each gene.}
#'     \item{p.values}{Numeric vector of corresponding p-values.}
#'   }
#'
#' @details
#' For each gene:
#' \itemize{
#'   \item Zero probabilities and dispersions are fixed from prior estimates.
#'   \item The model estimates group-specific mean parameters.
#'   \item When one group has all zero counts, a one-sided Z-test is applied instead.
#'   \item Otherwise, a standard two-sided Wald test is applied on the log-difference of group means.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' mat <- matrix(rnbinom(30, size = 5, mu = 5), nrow = 5)
#' grp <- c("A", "A", "A", "B", "B", "B")
#' y <- prepareDGE(mat, grp)
#' y <- sizeFactorsEst(y)
#' y <- tagwiseEst(y)
#' y <- hurdle_Wald_Test(y)
#' }
#'
#' @export
hurdle_Wald_Test <- function(y) {
  n_row <- nrow(y$counts)

  # Initialize output vectors
  wald_stats <- rep(NA, n_row)
  p_values   <- rep(1, n_row)

  groups_factor <- factor(y$samples$group)
  offsets       <- y$samples$size.factor
  levs          <- levels(groups_factor)

  for (i in 1:n_row) {
    current_counts <- as.numeric(y$counts[i, ])
    current_phi <- 1 / y$tagwise.disp[i]
    f_theta_A <- 1 - y$zero_prob_matrix[i, "prob_A"]
    f_theta_B <- 1 - y$zero_prob_matrix[i, "prob_B"]

    # Check for all-zero groups
    is_zero_A <- all(current_counts[groups_factor == levs[1]] <= 1)
    is_zero_B <- all(current_counts[groups_factor == levs[2]] <= 1)

    # Warm start
    m_A <- mean((current_counts[groups_factor == levs[1]] + 0.1) /
                  offsets[groups_factor == levs[1]], na.rm = TRUE)
    m_B <- mean((current_counts[groups_factor == levs[2]] + 0.1) /
                  offsets[groups_factor == levs[2]], na.rm = TRUE)

    tryCatch({
      fit <- optim(
        par = c(log(m_A), log(m_B) - log(m_A)),
        fn = nll_hurdle_fixed_theta_phi,
        counts = current_counts,
        groups = groups_factor,
        offset = offsets,
        fixed_theta_A = f_theta_A,
        fixed_theta_B = f_theta_B,
        fixed_phi = current_phi,
        one_group = FALSE,
        method = "BFGS",
        hessian = TRUE
      )

      H <- fit$hessian

      # Logic for One Group All Zeros
      if (is_zero_A || is_zero_B) {

        # choose which curvature to use
        idx <- if (is_zero_A) 2 else 1   # A all zero -> use delta direction; B all zero -> use A direction

        if (is.na(H[idx, idx]) || H[idx, idx] <= 1e-10) {
          wald_stats[i] <- 0
          p_values[i]   <- 1
        } else {
          se_beta <- sqrt(1 / H[idx, idx])

          # point estimate on log scale for the *non-zero group mean*
          beta_hat <- if (is_zero_A) (fit$par[1] + fit$par[2]) else fit$par[1]

          z_stat <- (beta_hat - 1) / se_beta
          wald_stats[i] <- z_stat
          p_values[i]   <- pnorm(z_stat, lower.tail = FALSE)  # one-sided
        }

      } else {
        # Logic for Standard Comparison (Neither is all zero)
        if (det(H) <= 1e-10 || any(is.na(H))) {
          wald_stats[i] <- 0
          p_values[i]   <- 1
        } else {
          vcov_mat <- solve(H)
          delta_beta <- fit$par[2]
          se_delta  <- sqrt(max(0, vcov_mat[2,2]))

          if (se_delta > 0) {
            z_stat <- delta_beta / se_delta
            wald_stats[i] <- z_stat
            # Two-sided for standard comparison
            p_values[i]   <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)
          } else {
            p_values[i] <- 1
          }
        }
      }

    }, error = function(e) {
      p_values[i] <- 1
    })
  }

  y$wald_stats <- wald_stats
  y$p.values   <- p_values
  return(y)
}

#' Interpolation Log-Likelihood with Diffeomorphisms
#'
#' Calculates log-likelihood and diffeomorphism map for a set of points
#'
#' @param X Numeric vector or nx1 matrix; represents data vector
#' @param Beta Numeric vector or mxp matrix; represents m different weight vectors for length p cosine basis; must have Euclidean norm less than pi.
#' @param b_vec Numeric vector; input points to interpolation; must be in order, first element equal to 0 and second element equal to 1.
#' @param lambda_vec Numeric vector; output points to interpolation; must have lambda_vec[i] > lambda_vec[i] < lambda_vec[i + 1]
#' @param interpolation 'cubic' or 'linear'; the interpolation method to use to fit b_vec & lambda_vec
#' @param transform_X boolean; whether to transform X into values in between 0 and 1 using min-max scaling; will do so automatically if min(X) < 0 or max(X) > 1
#' @param check_compat boolean; whether to perform process to check compatibility of given parameters
#' @returns List with following components: \item{diffeomorphism}{nxm matrix containing output after applying differomorphism to X for each beta weight vector} \item{log_like}{nxm matrix containing output after applying log-likelihood to X for each beta weight vector}
#' @export
#' @examples
#' # 3 modes (5 extrema) with two different beta weight vectors (-0.4, 0.3, 0.7) & (0.6, -0.1, 0.2)
#' X <- runif(100, 0, 1)
#' Beta <- matrix(c(-0.4, 0.3, 0.7, 0.6, -0.1, 0.2), nrow = 2, byrow = TRUE)
#' lambda_vec <- c(0, 1, 0.2, 0.4, 0.3, 2, 0)
#' b_vec <- seq(0, 1, length.out = length(lambda_vec))
#' result_lst <- diffeo_log_likelihood(X, Beta, b_vec, lambda_vec)
diffeo_log_likelihood <- function(X, Beta, b_vec, lambda_vec,
                                  interpolation = 'cubic', transform_X = TRUE,
                                  check_compat = TRUE) {

  # Turn X and Beta into matrices
  if (is.vector(X)){
    X <- as.matrix(X)
  }


  if (ncol(X) != 1){
    stop('X must be either a vector or a n x 1 matrix!')
  }

  # Transform X
  if ((transform_X) | (min(X) < 0) | (max(X) > 1)) {
    X <- min_max_transform_X(X)$new_X
  }


  # Check compatibility for interpolation and parameters
  Beta <- turn_1d_into_matrix(Beta)

  if ((interpolation != 'cubic') & (interpolation != 'linear')) {
    stop('interpolation must be either linear or cubic!')
  }

  if (check_compat) {
    Beta <- check_compat_log_lik(Beta, b_vec, lambda_vec)
  }

  # Apply gamma function to X and Beta and correct for floating point error
  apply_gamma <- gamma_function(X, Beta)
  apply_gamma[apply_gamma <= 0] <- 0
  apply_gamma[apply_gamma >= 1] <- 1

  # Apply interpolation to transformed space and transform back into matrix
  apply_gamma_vec <- as.vector(apply_gamma)
  interp_results_vec <- interp1(b_vec, lambda_vec, apply_gamma_vec,
                                method = interpolation)

  # This is to calculate normalizing constant
  # First getting interpolation results for a lot of points in (0, 1)
  interp_results_mat <- matrix(interp_results_vec, nrow = nrow(apply_gamma))
  apply_gamma_normalize <- gamma_function(X_FOR_INTEGRAL, Beta)
  apply_gamma_normalize[apply_gamma_normalize <= 0] <- 0
  apply_gamma_normalize[apply_gamma_normalize >= 1] <- 1

  apply_gamma_normalize_vec <- as.vector(apply_gamma_normalize)
  interp_normalize_vec <- interp1(b_vec, lambda_vec, apply_gamma_normalize_vec,
                                  method = interpolation)

  # Get normalizing constant by performing trapezoid integration
  interp_normalize_mat <- matrix(interp_normalize_vec,
                                 nrow = nrow(apply_gamma_normalize))
  normal_const <- apply(interp_normalize_mat, 2, TRAPZ_INTEGRAL_FUNC)

  # Get log likelihood
  log_like <- log(t(t(interp_results_mat)/normal_const))

  return(list(diffeomorphism = apply_gamma, log_like = log_like))

}

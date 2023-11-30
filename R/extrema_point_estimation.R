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
#' @returns list containing an nxm matrix containing output after applying differomorphism to X for each beta weight vector (diffeomorphism) and applying log-likelihood to X for each beta weight vector (log_like)
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

  if (abs(sum(Beta)) <= 1e-6){
    Beta <- Beta + matrix(runif(ncol(Beta), -0.001, 0.001), nrow = 1)
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


for_optimization_diffeo <- function(betas, input_X) {

  Beta <- matrix(betas, nrow = 1)

  log_like_full <- diffeo_log_likelihood(input_X, Beta, b_vec = c(0, 0.5, 1),
                                         lambda_vec = c(0, 1, 0),
                                         interpolation = 'cubic',
                                         transform_X = FALSE,
                                         check_compat = FALSE)
  return(-1 * sum(log_like_full$log_like))

}


hin1 <- function(x) {
  norms <- sum(x^2) - pi^2
  return(norms)
}

bayes_log_like <- function(betas, Y) {

  if (sum(betas^2) > pi^2) {
    return(-Inf)
  }

  Beta <- matrix(betas, nrow = 1)

  log_like_full <- diffeo_log_likelihood(input_X, Beta, b_vec = c(0, 0.5, 1),
                                         lambda_vec = c(0, 1, 0),
                                         interpolation = 'cubic',
                                         transform_X = FALSE,
                                         check_compat = FALSE)
  return(-1 * sum(log_like_full$log_like))
}

shit <- ess(bayes_log_like, input_X, Sig = 0.75 * diag(3), N_mcmc = 1000,
    burn_in = 1000, N = 3, FALSE)


plot(1:1000, shit[, 3], type = 'l')


global_optimize <- function(X, num_betas, num_trials = 25) {

  X <- rnorm(50)
  num_betas <- 3
  num_trials <- 25

  transformed_X <- min_max_transform_X(X)
  input_X <- transformed_X$new_X

  hist(input_X, breaks = 20)

  beta_bound <- pi/sqrt(num_betas)

  beta_starts <- matrix(runif(num_trials * num_betas, -beta_bound, beta_bound),
                        ncol = num_betas)

  fmin_value <- Inf

  for (j in 1:num_trials) {

    trial_point <- beta_starts[j, ]

    min_result <- fmincon(trial_point, for_optimization_diffeo, input_X = input_X,
                     method = 'SQP', A = NULL, b = NULL, Aeq = NULL,
                     beq = NULL, lb = NULL, ub = NULL, hin = hin1, heq = NULL,
                     tol = 1e-06, maxfeval = 10000, maxiter = 5000)

    if (min_result$value < fmin_value) {
      fmin_value <- min_result$value
      final_beta <- min_result$par
    }
  }

  X_test <- seq(0, 1, length.out = 2500)
  help <- diffeo_log_likelihood(X_test,
                                shit[950, ], b_vec = c(0, 0.5, 1),
                                lambda_vec = c(0, 1, 0),
                                interpolation = 'cubic',
                                transform_X = FALSE,
                                check_compat = FALSE)

  plot(X_test, exp(help$log_like), type = 'l')

  interp1(as.vector(help$diffeomorphism), X_test, 0.5)


}




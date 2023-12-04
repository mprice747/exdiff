#' Estimate Distribution and its Mode with Diffeomorphism MLE
#'
#' Obtain points estimates of pdf and mode by calculating MLE of diffeomorphism log-likelihood
#'
#' @param X Numeric vector; represents data vector
#' @param num_betas positive integer > 1; length of weight vector, also number of basis elements for cosine basis of diffeomorphism.
#' @param num_trials positive integer; number of random starting points, and thus trials, for optimization procedure.
#' @param beta_starts NULL or n x num_betas matrix; represents user chosen starting points. If NULL num_trials random beta vectors chosen.
#' @param plot_results boolean; whether to display plot which shows histogram of original data, the pdf estimate and the mode estimate.
#' @returns List with following components: \item{final_beta}{MLE weight estimate} \item{mode_estimate}{MLE mode estimate} \item{p_X}{input used for MLE pdf estimation} \item{mle_pdf}{MLE pdf estimate with p_X}
#' @export
#' @examples
#' # Sample from gamma(2, 1). Mode should be around 1
#' # X <- rgamma(50, 2, 1)
#' # mle_estimate <- diffeo_mle_estimate(X, num_betas = 3, num_trials = 25, plot_results = TRUE)
diffeo_mle_estimate <- function(X, num_betas, num_trials = 25, beta_starts = NULL,
                                plot_results = TRUE) {

  # num_betas need to be greater than or equal to 2
  check_integer_2(num_betas, 'num_betas')

  # num_trials need to be positive integer
  check_pos_integer(num_trials, 'num_trials')

  # Transform X to [0, 1] and obtain beta estimates
  transformed_X <- min_max_transform_X(X)
  input_X <- transformed_X$new_X

  final_beta <- global_optimize(input_X, num_betas, 'mle',
                                num_trials, beta_starts)

  # Uses estimates to calculate pdf
  fit_likelihood <- diffeo_log_likelihood(X_SAMPLE,
                                          final_beta, b_vec = c(0, 0.5, 1),
                                          lambda_vec = c(0, 1, 0),
                                          interpolation = 'cubic',
                                          transform_X = FALSE,
                                          check_compat = FALSE)

  # Get mode estimate by obtaining inverse of diffeomorphism
  mode_estimate_0_1 <- interp1(as.vector(fit_likelihood$diffeomorphism),
                               X_SAMPLE, 0.5)

  # Revert X back to original scale and recalculate mode and pdf
  range_X <- transformed_X$upper_bound - transformed_X$lower_bound
  mode_estimate <- mode_estimate_0_1 * range_X + transformed_X$lower_bound

  p_X <- X_SAMPLE * range_X + transformed_X$lower_bound
  pdf_unnormalized <- exp(fit_likelihood$log_like)

  mle_pdf <- as.vector(pdf_unnormalized/trapz(p_X, pdf_unnormalized))

  # Plot the MLE estimate of pdf and mode
  if (plot_results) {
    mle_plot <- mode_estimation_plot(X, p_X, mle_pdf, mode_estimate,
                                     'Mode Estimation Density Plot (MLE)')
    plot(mle_plot)
  }

  return(list(final_beta = final_beta, mode_estimate = mode_estimate,
              p_X = p_X, mle_pdf = pdf))
}

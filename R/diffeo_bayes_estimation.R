#' Estimate Distribution and Obtain Posterior of Mode with Diffeomorphism Bayesian Analysis
#'
#' Obtain MAP estimates of pdf and posterior distribution of modes by calculating posterior with adapative MCMC
#'
#' @param X Numeric vector; represents data vector.
#' @param num_betas positive integer > 1; length of weight vector, also number of basis elements for cosine basis of diffeomorphism.
#' @param num_trials positive integer; number of random starting points, and thus trials, for MAP optimization procedure.
#' @param beta_starts NULL or n x num_betas matrix; represents user chosen starting points for MAP. If NULL num_trials random beta vectors chosen.
#' @param num_samples positive integer; number of samples to return from MCMC.
#' @param burn_in positive integer; burn_in number for MCMC.
#' @param prior_sd positive real; prior is no covariance multivariate normal, represents standard deviation for each component.
#' @param prop_sd positive real; proposal distribution standard deviation for adaptive MCMC.
#' @param plot_results boolean; whether to display plot which shows histogram of original data, the pdf estimate and the mode estimate.
#' @returns List containing following components: \item{sampled_betas}{MCMC posterior sampled weight vectors} \item{sampled_modes}{MCMC posterior sampled modes} \item{bayes_map_beta}{MAP estimate of beta vector} \item{bayes_map_mode}{MAP estimate of mode} \item{p_X}{the input used for pdf estimation} \item{bayes_map_pdf}{MAP estimate of the pdf}
#' @export
#' @examples
#' # Sample from gamma(2, 1). Mode should be around 1
#' # X <- rgamma(50, 2, 1)
#' # bayes_estimate <- diffeo_bayes_estimate(X, num_betas = 3, num_samples = 5000, prior_sd = 0.75)
diffeo_bayes_estimate <- function(X, num_betas, num_trials = 25, beta_starts = NULL,
                                  num_samples = 5000, burn_in = 1000, prior_sd = 0.75,
                                   prop_sd = 0.5, plot_results = TRUE) {

  # num_betas need to be greater than or equal to 2
  check_integer_2(num_betas, 'num_betas')

  if (length(X) < num_betas){
    stop('num_betas must be less than the number of data points!')
  }

  # Following need to be positive integers
  check_pos_integer(num_trials, 'num_trials')
  check_pos_integer(num_samples, 'num_samples')
  check_pos_integer(burn_in, 'burn_in')

  # Following need to be positive numbers
  check_pos_number(prior_sd, 'prior_sd')
  check_pos_number(prop_sd, 'prop_sd')


  # Transform X to [0, 1] and obtain MAP estimates
  transformed_X <- min_max_transform_X(X)
  input_X <- transformed_X$new_X

  opt_type <- 'map'
  print('Starting Initial Guess Search (MAP)')
  map_estimate <- global_optimize(input_X, num_betas, opt_type, num_trials,
                                  beta_starts, prior_sd)
  print('Finished Initial Guess Search (MAP)')


  # Run MCMC sampler
  print('Starting Adaptive MCMC Sampler')
  bayes_estimates_0_1 <- mh_adaptive_sampling_betas(input_X, map_estimate, num_samples,
                                                    burn_in, prior_sd, prop_sd)
  print('Finished Adaptive MCMC Sampler')

  # Transform X back to original scale and obtain modified estimats
  range_X <- transformed_X$upper_bound - transformed_X$lower_bound

  sampled_betas <- bayes_estimates_0_1$sampled_betas
  sampled_modes <- bayes_estimates_0_1$sampled_modes * range_X + transformed_X$lower_bound
  bayes_map_mode <- bayes_estimates_0_1$bayes_map_mode * range_X + transformed_X$lower_bound
  p_X <- X_SAMPLE * range_X + transformed_X$lower_bound
  pdf_unnormalized <- exp(bayes_estimates_0_1$bayes_map_pdf)
  bayes_map_pdf <- as.vector(pdf_unnormalized/trapz(p_X, pdf_unnormalized))

  # Plot both the MAP estimates and the sampled modes
  if (plot_results) {
    bayes_plot <- diffeo_bayes_plot(X, p_X, bayes_map_pdf, bayes_map_mode, sampled_modes)
    plot(bayes_plot)
  }

  return(list(sampled_betas = sampled_betas, sampled_modes = sampled_modes,
              bayes_map_beta = map_estimate,
              bayes_map_mode = bayes_map_mode, p_X = p_X,
              bayes_map_pdf = bayes_map_pdf))
}

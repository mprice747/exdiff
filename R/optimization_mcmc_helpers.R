# Helpers for MLE optimization and MCMC sampler

#' Constant used for plotting likelihood functions anf getting mode estimates
#' @noRd
X_SAMPLE <- seq(0, 1, length.out = 1000)

#' Given beta (weight) vector and dataset with values in between 0 and 1, calculate negative joint log likelihood (Only one mode).
#' @param betas numeric vector; represents weight vector for cosine basis
#' @param input_X  numeric vector; represents data for log likelihood. Points must be in between 0 and 1
#' @returns negative log likelihood
#' @noRd
mle_optimization <- function(betas, input_X) {


  Beta <- matrix(betas, nrow = 1)

  # If norm > pi, return infinity (invalid point)
  if (sum(Beta^2) >(pi^2) + 1e-05) {
    return(Inf)
  }

  # Get negative joint log-likelihood given beta and input_X
  log_like_full <- diffeo_log_likelihood(input_X, Beta, b_vec = c(0, 0.5, 1),
                                         lambda_vec = c(0, 1, 0),
                                         interpolation = 'cubic',
                                         transform_X = FALSE,
                                         check_compat = FALSE)

  return(-1 * sum(log_like_full$log_like))
}


#' Function used for fmincon. Makes sure norm of betas < pi
#' @param x vector
#' @returns difference of Euclidean norm^2 and pi^2
#' @noRd
hin1 <- function(x) {
  norms <- sum(x^2) - pi^2
  return(norms)
}

#' Check if positive number
#' @param x number to check
#' @param name variable name of number to check
#' @noRd
check_pos_number <- function(x, name) {
  if (x <= 0){
    stop(paste(name, 'needs to be a positive number!'))
  }
}

#' Check if positive integer
#' @param x number to check
#' @param name variable name of number to check
#' @noRd
check_pos_integer <- function(x, name){
  if ((x < 1) | (x%%1 != 0) ){
    stop(paste(name, 'needs to be a positive integer!'))
  }
}

#' Check if positive integer greater than or equal to 2
#' @param x number to check
#' @param name variable name of number to check
#' @noRd
check_integer_2 <- function(x, name){
  if ((x < 2) | (x%%1 != 0) ){
    stop(paste(name, 'needs to be an integer greater than 2!'))
  }
}
#' Given dataset, find optimal beta vector that maximizes MLE or MAP (Only assumes one mode, so designed to find most extreme mode in probability distribution) Uses fmincon with multiple starting points and picks
#' @param input_X numeric vector; represents data for log likelihood. Points must be in between 0 and 1
#' @param num_betas positive integer; length of weight vector, also number of basis elements for cosine basis of diffeomorphism
#' @param optimize_type either 'mle' or 'map', whether to find mle or map estimate
#' @param num_trials positive integer; number of random starting points for fmincon
#' @param beta_starts - NULL or n x num_betas matrix; represents user chosen starting points. If NULL num_trials random beta vectors chosen
#' @param prior_sd positive real or NULL; Represents prior distribution standard deviation for MAP
#' @returns numeric vector representing fmincon's estimate of the argmin should be
#' @noRd
global_optimize <- function(input_X, num_betas, optimize_type,
                            num_trials = 25, beta_starts = NULL, prior_sd = NULL) {

  # If beta_starts is null generate random starting points
  if (is.null(beta_starts)) {

    beta_bound <- pi/sqrt(num_betas)
    beta_starts <- matrix(runif(num_trials * num_betas, -beta_bound, beta_bound),
                          ncol = num_betas)
  } else{
    # If beta_starts is provided, checks if matches with num_betas
    if (num_betas != ncol(beta_starts)) {
      stop('num_betas does not match the number of columns in beta_starts!')
    }

  }

  if (optimize_type == 'map') {

    if(is.null(prior_sd)) {
      stop('prior_sd can not be null if using map estimate!')
    }

    # If MAP, initializes cost function used for minimization
    prior_sd_vec <- rep(prior_sd, num_betas)
    prior_mean_vec <- rep(0, num_betas)

    map_optimization <- function(betas, input_X) {
      return( mle_optimization(betas, input_X) - sum(dnorm(betas,
                  mean = prior_mean_vec, sd = prior_sd_vec, log = TRUE) ))
    }

    func_to_minimize <- map_optimization

  } else if (optimize_type == 'mle') {

    # If MLE, uses mle_optimization for cost function
    func_to_minimize <- mle_optimization

  } else{
    stop('optimize_type must be mle or map!')
  }

  fmin_value <- Inf

  # For each trial point, run fmincon and see which beta minimized the cost function
  # the most
  for (j in 1:num_trials) {

    trial_point <- beta_starts[j, ]

    # Possible fmincon errors, just skip provided particular beta_start
    tryCatch(expr = {min_result <- fmincon(trial_point,
                                           func_to_minimize, input_X = input_X,
                                           method = 'SQP', A = NULL, b = NULL, Aeq = NULL,
                                           beq = NULL, lb = NULL, ub = NULL, hin = hin1, heq = NULL,
                                           tol = 1e-06, maxfeval = 10000, maxiter = 5000)
    },error=function(e) {
      message(paste('Invalid beta vector for fmincon, skipping ', j, 'th beta_start vector',
                    sep = ''))
    })


    if (min_result$value < fmin_value) {
      fmin_value <- min_result$value
      final_beta <- min_result$par
    }

  }

  # Return best result
  return(final_beta)

}

#' For Bayesian estimation, runs adaptive MCMC to get posterior distribution of beta_vectors
#' @param input_X numeric vector; represents data for log likelihood. Points must be in between 0 and 1
#' @param map_estimate NULL numeric vector; MAP estimate which is starting point for MCMC
#' @param num_samples positive integer; number of samples to return from MCMC
#' @param burn_in positive integer; burn_in number for MCMC
#' @param prior_sd positive real; prior is no covariance multivariate normal, represents standard deviation for each component
#' @param prop_sd positive real; proposal distribution standard deviation for adaptive MCMC
#' @returns list containing sampled beta vectors (sampled_betas), sampled modes calculated from the sampled beta vectors (sampled_modes), MAP estimate of the mode (bayes_map_mode), and MAP estimate of the pdf using X_SAMPLE as input (bayes_map_pdf)
#' @noRd
mh_adaptive_sampling_betas <- function(input_X, map_estimate, num_samples = 5000,
                                       burn_in = 1000, prior_sd = 0.75,
                                       prop_sd = 0.5) {

  # Intialize vectors for MCMC function
  p <- length(map_estimate)
  prior_sd_vec <- rep(prior_sd, p)
  prior_mean_vec <- rep(0, p)

  function_for_mcmc <- function(betas) {
    return(-1 * mle_optimization(betas, input_X) + sum(dnorm(betas,
         mean = prior_mean_vec, sd = prior_sd_vec, log = TRUE) ))
  }

  # Run Adaptive MCMC and get posterior results
  initial_scale <- rep(prop_sd, p)
  adaptive_MCMC <- MCMC(p = function_for_mcmc, init = map_estimate,
                        n = num_samples + burn_in, scale = initial_scale, adapt = TRUE,
                        acc.rate = 0.3)

  # Get sampled betas
  sampled_betas <- adaptive_MCMC$samples[(burn_in + 1):(num_samples + burn_in), ]

  # Get all likelihoods to calculate estimate mode (inverse of diffeomorphism applied at 0.5)
  all_likelihoods <- diffeo_log_likelihood(X_SAMPLE, sampled_betas, b_vec = c(0, 0.5, 1),
                                           lambda_vec = c(0, 1, 0),
                                           interpolation = 'cubic',
                                           transform_X = FALSE,
                                           check_compat = FALSE)

  sampled_modes <- apply(all_likelihoods$diffeomorphism, 2,
                         interp1, y = X_SAMPLE, xi = 0.5)

  # Get map mode using map estimate
  bayes_map_results <- diffeo_log_likelihood(X_SAMPLE, map_estimate,
                                             b_vec = c(0, 0.5, 1),
                                             lambda_vec = c(0, 1, 0),
                                             interpolation = 'cubic',
                                             transform_X = FALSE,
                                             check_compat = FALSE)

  bayes_map_mode <- interp1(as.vector(bayes_map_results$diffeomorphism),
                            X_SAMPLE, 0.5)

  return(list(sampled_betas = sampled_betas, sampled_modes = sampled_modes,
              bayes_map_mode = bayes_map_mode,
              bayes_map_pdf = as.vector(exp(bayes_map_results$log_like))))
}

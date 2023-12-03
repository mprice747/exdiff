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

X_SAMPLE <- seq(0, 1, length.out = 1000)

mle_optimization <- function(betas, input_X) {

  Beta <- matrix(betas, nrow = 1)

  if (sum(Beta^2) >(pi^2) + 1e-05) {
    return(Inf)
  }

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

global_optimize <- function(input_X, num_betas, optimize_type,
                            num_trials = 25, beta_starts = NULL, prior_sd = NULL) {

  if (is.null(beta_starts)) {
    beta_bound <- pi/sqrt(num_betas)
    beta_starts <- matrix(runif(num_trials * num_betas, -beta_bound, beta_bound),
                          ncol = num_betas)
  }

  if(optimize_type == 'map') {
    prior_sd_vec <- rep(prior_sd, num_betas)
    prior_mean_vec <- rep(0, num_betas)

    map_optimization <- function(betas) {
      return( mle_optimization(betas, input_X) - sum(dnorm(betas,
            mean = prior_mean_vec, sd = prior_sd_vec, log = TRUE) ))
    }
    func_to_minimize <- map_optimization
  }
  else{
    func_to_minimize <- mle_optimization
  }


  fmin_value <- Inf

  for (j in 1:num_trials) {

    trial_point <- beta_starts[j, ]


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

  return(final_beta)

}

mode_estimation_plot <- function(X, p_X, pdf, mode_estimate) {

  where_mode <- mode_estimate/max(p_X)

  if (where_mode <= 0.5){
    annotate_x <- Inf
    annotate_hjust <- 1
  }
  else {
    annotate_x <- -Inf
    annotate_hjust <- 0
  }

  mode_pdf <- interp1(p_X, pdf, mode_estimate)

  plot_mle <- ggplot() + geom_histogram(aes(X, after_stat(density)), colour = 1, fill = 'white',
              bins = as.integer(length(X) * 1/2.5) ) + theme(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank()) + geom_line(aes(p_X, pdf),
            size = 1.5) + geom_segment(aes(x = mode_estimate, y = 0, xend = mode_estimate,
            yend = mode_pdf), color = "red", size=1.5) + annotate("label", x = annotate_x, y = Inf, label = paste("Mode =", round(mode_estimate, 3)),
           hjust = annotate_hjust, vjust = 1, size = 5.5) + xlab('Input Data') + ylab('Density Estimation') + ggtitle('Mode Estimation Density Plot')

  return(plot_mle)
}

diffeo_mle_estimate <- function(X, num_betas, num_trials = 25, beta_starts = NULL,
                                plot_results = TRUE) {

  transformed_X <- min_max_transform_X(X)
  input_X <- transformed_X$new_X

  final_beta <- global_optimize(input_X, num_betas, 'mle', num_trials, beta_starts)

  fit_likelihood <- diffeo_log_likelihood(X_SAMPLE,
                                final_beta, b_vec = c(0, 0.5, 1),
                                lambda_vec = c(0, 1, 0),
                                interpolation = 'cubic',
                                transform_X = FALSE,
                                check_compat = FALSE)

  mode_estimate_0_1 <- interp1(as.vector(fit_likelihood$diffeomorphism),
                              X_SAMPLE, 0.5)

  range_X <- transformed_X$upper_bound - transformed_X$lower_bound

  mode_estimate <- mode_estimate_0_1 * range_X + transformed_X$lower_bound

  p_X <- X_SAMPLE * range_X + transformed_X$lower_bound
  pdf_unnormalized <- exp(fit_likelihood$log_like)

  pdf <- as.vector(pdf_unnormalized/trapz(p_X, pdf_unnormalized))

  if (plot_results) {
    mle_plot <- mode_estimation_plot(X, p_X, pdf, mode_estimate)
    plot(mle_plot)
  }

  return(list(final_beta = final_beta, mode_estimate = mode_estimate,
              p_X = p_X, pdf = pdf ))
}


mh_adaptive_sampling_betas <- function(input_X, initial_guess, num_samples = 5000,
                                       burn_in = 1000, prior_sd = 0.75,
                                       prop_sd = 0.5) {

  p <- length(initial_guess)

  prior_sd_vec <- rep(prior_sd, p)
  prior_mean_vec <- rep(0, p)

  function_for_mcmc <- function(betas) {
    return(-1 * mle_optimization(betas, input_X) + sum(dnorm(betas,
                                                             mean = prior_mean_vec, sd = prior_sd_vec, log = TRUE) ))
  }


  initial_scale <- rep(prop_sd, p)

  adaptive_MCMC <- MCMC(p = function_for_mcmc, init = initial_guess,
                        n = num_samples + burn_in, scale = initial_scale, adapt = TRUE,
                        acc.rate = 0.3)

  sampled_betas <- adaptive_MCMC$samples[(burn_in + 1):(num_samples + burn_in), ]
  bayes_avg_betas <- colMeans(sampled_betas)

  all_likelihoods <- diffeo_log_likelihood(X_SAMPLE, sampled_betas, b_vec = c(0, 0.5, 1),
                                           lambda_vec = c(0, 1, 0),
                                           interpolation = 'cubic',
                                           transform_X = FALSE,
                                           check_compat = FALSE)

  sampled_modes <- apply(all_likelihoods$diffeomorphism, 2,
                             interp1, y = X_SAMPLE, xi = 0.5)


  bayes_avg_results <- diffeo_log_likelihood(X_SAMPLE, bayes_avg_betas,
                                                b_vec = c(0, 0.5, 1),
                                                lambda_vec = c(0, 1, 0),
                                                interpolation = 'cubic',
                                                transform_X = FALSE,
                                                check_compat = FALSE)

  bayes_avg_mode <- interp1(as.vector(bayes_avg_results$diffeomorphism),
                               X_SAMPLE, 0.5)

  return(list(sampled_betas = sampled_betas, sampled_modes = sampled_modes,
              bayes_avg_mode = bayes_avg_mode,
              bayes_avg_likelihood = as.vector(exp(bayes_avg_results$log_like))))
}

diffeo_bayes_plot <- function(X, p_X, pdf, mode_estimate, sampled_modes) {

  plot_1 <- mode_estimation_plot(X, p_X, pdf, mode_estimate)

  plot_2 <- ggplot() + geom_histogram(aes(sampled_modes, after_stat(density)), colour = 1, fill = 'white',
           bins =20) + theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + xlab('Sampled Modes') + ylab('Density') + ggtitle('Sampled Modes Histogram')

  return(grid.arrange(plot_1, plot_2, ncol=2))
}

diffeo_bayes_estimates <- function(X, num_betas, num_samples = 5000,
                                   burn_in = 1000, prior_sd = 0.75,
                                   prop_sd = 0.5, plot_results = TRUE) {

  transformed_X <- min_max_transform_X(X)
  input_X <- transformed_X$new_X


  print('Starting Initial Guess Search (MAP)')
  initial_guess <- global_optimize(input_X, num_betas, 'map', 25, NULL, prior_sd)
  print('Finished Initial Guess Search (MAP)')


  print('Starting Adaptive MCMC Sampler')
  bayes_estimates_0_1 <- mh_adaptive_sampling_betas(input_X, initial_guess,num_samples,
                             burn_in, prior_sd, prop_sd)
  print('Finished Adaptive MCMC Sampler')

  range_X <- transformed_X$upper_bound - transformed_X$lower_bound

  sampled_betas <- bayes_estimates_0_1$sampled_betas
  sampled_modes <- bayes_estimates_0_1$sampled_modes * range_X + transformed_X$lower_bound
  bayes_avg_mode <- bayes_estimates_0_1$bayes_avg_mode * range_X + transformed_X$lower_bound

  p_X <- X_SAMPLE * range_X + transformed_X$lower_bound
  pdf_unnormalized <- exp(bayes_estimates_0_1$bayes_avg_likelihood)
  bayes_pdf <- as.vector(pdf_unnormalized/trapz(p_X, pdf_unnormalized))

  if (plot_results) {
    bayes_plot <- diffeo_bayes_plot(X, p_X, bayes_pdf, bayes_avg_mode, sampled_modes)
    plot(bayes_plot)
  }

  return(list(sampled_betas = sampled_betas, sampled_modes = sampled_modes,
              bayes_avg_mode = bayes_avg_mode, p_X = p_X,
              bayes_pdf = bayes_pdf))
}




set.seed(42)
X <- rgamma(50, 2, 1)
hist(X, breaks = 10)

input_X <- min_max_transform_X(X)$new_X

hist(input_X)

initial_guess <- global_optimize(input_X, 3, 25)
mle_estimate <- diffeo_mle_estimate(X, 3)


bayes_estimators <- diffeo_bayes_estimates(X, 3)

mean(bayes_estimators$bayes_avg_mode)


plot(bayes_estimators$p_X, bayes_estimators$bayes_pdf, type = 'l')










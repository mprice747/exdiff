# Helpers for plotting functionalities

#' Displays histogram of original data, pdf estimate and mode estimate
#' @param X numeric vector; original data
#' @param p_X numeric vector; data used for pdf estimation
#' @param pdf numeric vector; pdf output with p_X as input
#' @param mode_estimate real number; estimated mode, to be displayed on top left or top right corner
#' @param plot_title string; title of plot
mode_estimation_plot <- function(X, p_X, pdf, mode_estimate, plot_title) {

  # Whether mode should be displayed on top left or top right
  where_mode <- mode_estimate/max(p_X)
  if (where_mode <= 0.5){
    annotate_x <- Inf
    annotate_hjust <- 1
  } else {
    annotate_x <- -Inf
    annotate_hjust <- 0
  }

  # Get pdf estimate of the mode
  mode_pdf <- interp1(p_X, pdf, mode_estimate)

  # Make and return plot
  estimate_plot <- ggplot() + geom_histogram(aes(X, after_stat(density)), colour = 1, fill = 'white',
                    bins = as.integer(length(X) * 1/2.5) ) + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()) + geom_line(aes(p_X, pdf),
                     size = 1.5) + geom_segment(aes(x = mode_estimate, y = 0, xend = mode_estimate,
                   yend = mode_pdf), color = "red", size=1.5) + annotate("label", x = annotate_x, y = Inf, label = paste("Mode =", round(mode_estimate, 3)),
      hjust = annotate_hjust, vjust = 1, size = 5.5) + xlab('Input Data') + ylab('Density Estimation') + ggtitle(plot_title)

  return(estimate_plot)
}

#' Plot for bayesian analysis of modes. Displays mode estimation plot with MAP estimates on top and histogram of MCMC sampled modes on bottom
#' @param X numeric vector; original data
#' @param p_X numeric vector; data used for pdf estimation
#' @param pdf numeric vector; pdf output with p_X as input
#' @param mode_estimate real number; estimated mode, to be displayed on top left or top right corner
#' @param sampled_modes numeric vector; contains sampled modes from Bayesian MCMC
diffeo_bayes_plot <- function(X, p_X, pdf, mode_estimate, sampled_modes) {

  #

  # Mode estimation plot
  plot_1 <- mode_estimation_plot(X, p_X, pdf, mode_estimate, 'Mode Estimation Density (MAP)')

  # Histogram of sampled modes
  plot_2 <- ggplot() + geom_histogram(aes(sampled_modes, after_stat(density)), colour = 1, fill = 'white',
                          bins =20) + theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank()) + xlab('Sampled Modes') + ylab('Density') + ggtitle('Sampled Modes Histogram')

  return(grid.arrange(plot_1, plot_2, nrow =2))
}


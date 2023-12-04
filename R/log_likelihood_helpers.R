# Helpers to calculate diffeomorphism log_likelihood

#' SQRT(2) constant
#' @noRd
SQRT_2 <- sqrt(2)

#' Used for calculating normalizing constant (Trapezoid integration)
#' @noRd
X_FOR_INTEGRAL <- as.matrix(seq(0, 1, length.out = 2000))

#' Calculating integral from 0 to 1 of function via trapezoid integral
#' @param func function to integrate
#' @returns integral result
#' @noRd
TRAPZ_INTEGRAL_FUNC <- function(func) {
  return(trapz(X_FOR_INTEGRAL,func))
}


#' Check if input is matrix. If not, turn into 1 column matrix
#' @param X array like
#' @returns matrix form of X (nx1)
#' @noRd
check_matrix <- function(X) {

  if (!is.matrix(X)) {
    X <- matrix(X, ncol = 1)
  }

  return(X)
}

#' Check if input is matrix. If not, turn into 1 row matrix
#' @param X array like
#' @returns matrix form of X (1xn)
#' @noRd
turn_1d_into_matrix <- function(X) {

  if (is.null(dim(X))) {
    X <- matrix(X, nrow = 1)
    return(X)
  }
  return(X)
}

#' Given weights (Beta), calculates cosine basis evaluated at X
#' @param X nx1 Matrix of observations
#' @param Beta mxj Matrix of different beta vectors (weight parameters)
#' @returns nxm Matrix of representing the n observations evaluated at m different cosine bases
#' @noRd
cosine_basis <- function(X, Beta) {

  Cos_Mat <- cos(pi * (X %*% matrix(seq(1, ncol(Beta)), nrow = 1)))

  result <- tcrossprod(Cos_Mat, (SQRT_2 * Beta))

  return(result)
}

#' Calculate integral from 0 to X of cosine basis with beta as weights
#' @param X nx1 Matrix of observations description
#' @param Beta mxj Matrix of different beta vectors (weight parameters)
#' @returns nxm Matrix of representing the integral of 0 to X of cosine basis with different weight parameters
#' @noRd
int_cosine_basis <- function(X, Beta) {

  j_vec <- seq(1, ncol(Beta))

  Sin_Mat <- sin(pi * (X %*% matrix(j_vec, nrow = 1)))

  result <- Sin_Mat %*% (t((SQRT_2 * Beta))/(j_vec * pi))

  return(result)
}

#' Calculate int_{0}^{x}cos^2(y * pi * j)dy
#' @param x number, upper bound
#' @param j positive integer
#' @returns integral result
#' @noRd
cos_2_integral <- function(x, j) {

  two_pi_j <- 2 * pi * j
  result <- sin(two_pi_j * x)/(2 * two_pi_j) + x/2

  return(result)

}
#' Calculate int_{0}^{x} cos(y * pi * j1)cos(y * pi * j2)dy where j1 != j2
#' @param x number, upper bound
#' @param j1 positive integer
#' @param j2 positive integer, different from j1
#' @returns integral result
#' @noRd
cos_cos_integral <- function(x, j1, j2){

  pi_j1_j2 <- pi * (j1 + j2)
  pi_j2_sub_j1 <- pi * (j2 - j1)

  result <- sin(pi_j1_j2 * x)/(2 * pi_j1_j2) +
    sin(pi_j2_sub_j1 * x)/(2 * pi_j2_sub_j1)

  return(result)
}

#' Given vector, output 2 x choose(length(vec), 2) containing every combination of elements within vec
#' @param vec 1-d vector
#' @returns vec containing every possible product of two numbers
#' @noRd
different_beta_product <- function(vec) {

  vec_comb <- combn(vec, 2)

  return(vec_comb[1, ] * vec_comb[2, ])
}

#' Calculate integral from 0 to X of cosine basis squared with beta as weights
#' @param X nx1 Matrix of observations
#' @param Beta mxj Matrix of different beta vectors (weight parameters)
#' @returns nxm Matrix of representing the integral of 0 to X of cosine basis squared with different weight parameters
#' @noRd
int_cosine_basis_2 <- function(X, Beta) {

  # Produce every possible pair of position
  j_seq <- seq(1, ncol(Beta))
  j_comb <- combn(j_seq, 2)


  # Get cos^2 squared integral of different basis functions
  cos_cos_ints <- turn_1d_into_matrix(apply(X, 1, cos_cos_integral,
                                            j1 = j_comb[1, ], j2 = j_comb[2, ]))
  diff_betas_mult <- turn_1d_into_matrix(apply(Beta, 1, different_beta_product))
  cos_cos_all <- 4 * crossprod(cos_cos_ints, diff_betas_mult)

  # Get cos^2 squared integral of same basis functions
  cos_2_ints <- turn_1d_into_matrix(apply(X, 1, cos_2_integral, j = j_seq))

  # Add all to return result
  cos_2_all <- 2 * t(Beta^2 %*% cos_2_ints)
  result <- cos_cos_all + cos_2_all

  return(result)
}

#' Gamma function which transforms domain (0, 1) to new domain (0, 1) via diffeomorphism
#' @param X nx1 Matrix of observations
#' @param Beta  mxj Matrix of different beta vectors (weight parameters)
#' @returns nxm Matrix of representing the results of the gamma function defined in the paper (int_{0}^{X} exponential_map(t)^2 dt)
#' @noRd
gamma_function <- function(X, Beta) {


  # Find integral of cosine basis and cosine basis squared
  Beta_norms <- sqrt(rowSums(Beta^2))

  int_1 <- t(int_cosine_basis(X, Beta))

  int_2 <- t(int_cosine_basis_2(X, Beta))

  # cosine and sines of norms in exponential map
  cos_beta <- cos(Beta_norms)
  sin_beta_div <- sin(Beta_norms)/Beta_norms
  cos_beta_X <- outer(cos(Beta_norms)^2, X[, 1])

  # Full integral of square exponential map
  exp_map <- cos_beta_X + ((2 * cos_beta *
                              sin_beta_div) * int_1) + sin_beta_div^2 * int_2

  return(t(exp_map))

}


#' Helper to transform X to interval in between 0 and 1 (Required for diffeomorphism)
#' @param X vector or matrix containing data
#' @returns list containing transformed X (new_X), lower_bound used in transformation and upper bound used in transformation
#' @noRd
min_max_transform_X <- function(X) {


  sd_X <- sd(X)
  sqrt_n <- sqrt(length(X))

  # Bounds are at extreme points plus or minus standard error
  lower_bound <- min(X) - sd_X/sqrt_n
  upper_bound <- max(X) + sd_X/sqrt_n

  # Min-max scaling
  new_X <- (X - lower_bound)/(upper_bound - lower_bound)

  return(list(new_X = new_X, lower_bound = lower_bound, upper_bound = upper_bound))

}

#' Check if in lambda vector lambda_{i - 1} > lambda{i} < lambda_{i + 1} when i is odd (lambda_vec length 5 or higher)
#' @param lambda_vec vector of positive numbers
#' @noRd
check_fluctuating_lambda <- function(lambda_vec) {

  check_fluctuating <- seq(from = 3, to = length(lambda_vec) - 2, by = 2)
  for (valley in check_fluctuating) {
    if ((lambda_vec[valley] >= lambda_vec[valley - 1]) | (lambda_vec[valley] >= lambda_vec[valley + 1])  ){
      stop('At least one lambda in lambda_vec misspecified! (not lambda_{i - 1} > lambda{i} < lambda_{i + 1}) ')
    }
  }
}

#' # Check compatibility of all parameters of log-likelihood function
#' @param Beta Numeric vector or mxp matrix; represents m different weight vectors for length p cosine basis; must have Euclidean norm less than pi.
#' @param b_vec Numeric vector; input points to interpolation; must be in order, first element equal to 0 and second element equal to 1.
#' @param lambda_vec Numeric vector; output points to interpolation; must have lambda_vec[i] > lambda_vec[i] < lambda_vec[i + 1].
#' @returns modified beta matrix in case some needed to be removed
#' @noRd
check_compat_log_lik <- function(Beta, b_vec, lambda_vec){



  # Both b_vec and lambda_vec must have the same, odd length
  if (length(b_vec) != length(lambda_vec)) {
    stop('b_vec and lambda_vec do not have the same length!')
  }

  if ((length(b_vec) %% 2) == 0 ){
    stop('b_vec and lambda_vec must have odd number length!')
  }

  # b_vec must be sorted and have its first element equal to 0 and last equal to 1
  if ((b_vec[1] != 0) | (b_vec[length(b_vec)] != 1)){
    stop('b_vec\'s first and last elements must equal 0 and 1 respectively!')
  }
  if (is.unsorted(b_vec)){
    stop('b_vec is unsorted!')
  }

  # lambda_vec's first and last elements equal to 0
  if ((lambda_vec[1] != 0) | (lambda_vec[length(lambda_vec)] != 0)) {
    stop('lambda_vec\'s first and last elements must equal 0')
  }

  # lambda_vec's second element must equal 1
  if ((lambda_vec[2] != 1)) {
    stop('lambda_vec\'s second element must equal 1')
  }

  # Check if lambda_vec is "fluctuating" if length greater than 3
  if (length(lambda_vec) > 3) {
    check_fluctuating_lambda(lambda_vec)
  }

  if(ncol(Beta) == 1){
    stop('Beta vectors must have bigger length than 1!')
  }

  # Check if any beta norm is greater than pi^2. If so, remove from calculation
  beta_norms <- rowSums(Beta^2)
  incompat_beta <- which((beta_norms >= pi^2) | (beta_norms <= 1e-6))

  if (length(incompat_beta) == length(beta_norms)) {
    stop('Norms of all supplied weight vectors > pi or near 0! Stopping log likelihood calculation')
  }
  if (length(incompat_beta) > 0){
    print('Subset of weight vectors have norm > pi or near 0! Removing them from log likelihood calculation')
    new_Beta <- Beta[-incompat_beta, ]
    return(new_Beta)
  }

  return(Beta)
}



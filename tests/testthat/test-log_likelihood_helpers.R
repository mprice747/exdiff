# Make Sure Diffeomorphism produces transform from [0, 1] to [0, 1] and still increasing
X <- as.matrix(seq(0, 1, length.out = 100))
Beta_1 <- matrix(c(0.5, -0.5, 0.43), nrow = 1)
Beta_2 <- matrix(c(1, -0.6), nrow = 1)
Beta_3 <- matrix(c(0.01, -0.01, 0.02, -0.02), nrow = 1)

# Function to check if vector or 1d matrix is increasing
is_increasing <- function(test_vec) {
  n <- length(test_vec)
  if (n == 1){
    return(T)
  }

  num_increase <- sum(test_vec[2:n] - test_vec[1:(n - 1)] >= 0)
  return(num_increase == (n - 1))
}

# Test if gamma function produces gamma(0) = 0, gamma(1) = 1 and produces increasing sequence
# X MUST BE INCREASING
test_gamma_function <- function(X, Beta) {

  # Get gamma function result
  n <- length(X)
  gamma_mat <- gamma_function(X, Beta)

  # Gamma output increasing
  gamma_increasing <- is_increasing(gamma_mat)

  gamma_list <- list(first_ele = gamma_mat[1], last_ele = gamma_mat[n],
                     gamma_increasing = gamma_increasing )

  return(gamma_list)
}

gamma_list_1 <- test_gamma_function(X, Beta_1)
gamma_list_2 <- test_gamma_function(X, Beta_2)
gamma_list_3 <- test_gamma_function(X, Beta_3)

# Gamma(0) == 0
test_that('Test if gamma(0) == 0',
          {expect_lte(abs(gamma_list_1$first_ele), 1e-6)
            expect_lte(abs(gamma_list_2$first_ele), 1e-6)
            expect_lte(abs(gamma_list_3$first_ele), 1e-6)})

# Gamma(1) == 1
test_that('Test if gamma(1) == 1',
          {expect_lte(abs(gamma_list_1$last_ele - 1), 1e-6)
            expect_lte(abs(gamma_list_2$last_ele - 1), 1e-6)
            expect_lte(abs(gamma_list_3$last_ele - 1), 1e-6)})

# Increasing gamma values
test_that('Test if gamma function is producing increasing values from 0 to 1',
          {expect_true(gamma_list_1$gamma_increasing)
            expect_true(gamma_list_2$gamma_increasing)
            expect_true(gamma_list_3$gamma_increasing)})

set.seed(42)
X_transform_1 <- rnorm(50)
X_transform_2 <- runif(100, -500, -200)
X_transform_3 <- rgamma(1000, 1)

# Test if min_max_transform_X transform X to be in between 0 and 1
test_min_max_transform_X <- function(X_transform){
  new_X <- min_max_transform_X(X_transform)$new_X
  return((min(new_X) > 0) & (max(new_X) < 1))
}

test_that('Test if min_max_transform_X transforms data to be in between 0 and 1',
          {expect_true(test_min_max_transform_X(X_transform_1))
            expect_true(test_min_max_transform_X(X_transform_2))
            expect_true(test_min_max_transform_X(X_transform_3))})


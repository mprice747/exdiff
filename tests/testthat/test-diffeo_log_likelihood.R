# Test the diffeo log likelihood function
X <- as.matrix(seq(0.01, 0.99, length.out = 2000))

# Initialize 3 sets of parameters (1, 3 and 5 extrema)
Beta_1 <- matrix(c(0.41, -0.99, 0.76), nrow = 1)
Beta_2 <- matrix(c(1, -1.5), nrow = 1)
Beta_3 <- matrix(c(0.05, -0.05, 0.02, -0.01), nrow = 1)

b_vec_1 <- seq(0, 1, length.out = 3)
b_vec_2 <- seq(0, 1, length.out = 5)
b_vec_3 <- seq(0, 1, length.out = 7)

lambda_vec_1 <- c(0, 1, 0)
lambda_vec_2 <- c(0, 1, 0.5, 0.76, 0)
lambda_vec_3 <- c(0, 1, 0.75, 2, 0.54, 1.25, 0)

# Given a vector find number of extrema by looping through elements
find_number_extrema <- function(Y) {

  num_extrema <- 0
  extrema_location <- c()

  # See if Y[i] is a peak or valley
  for (i in 2:(length(Y) - 1)) {

    peak <- (Y[i] > Y[i - 1]) & (Y[i] > Y[i + 1])
    valley <- (Y[i] < Y[i - 1]) & (Y[i] < Y[i + 1])

    # If so add to list
    if (peak | valley) {
      num_extrema <- num_extrema + 1
      extrema_location <- append(extrema_location, i)
    }
  }
  # Return number of extrema and location
  return(list(num_extrema = num_extrema, extrema_location = extrema_location))
}

# Test diffeomorphism log likelihood function
# Check if number of extrema remain the same as inputted, the height ratios are the same
# inputted and density area is around 1
test_diffeo_log_like <- function(X, Beta, b_vec, lambda_vec){

  diffeo <- diffeo_log_likelihood(X, Beta, b_vec, lambda_vec,
                                    transform_X = FALSE)

  log_like <- as.vector(diffeo$log_like)

  # Find extrema and their location
  found_extrema <- find_number_extrema(log_like)
  extrema_location <- found_extrema$extrema_location
  extrema_values <- exp(log_like[extrema_location])

  # Get error in lambda vector
  llv_1 <- length(lambda_vec) - 1
  extrema_heights_diff <- mean(abs((extrema_values/extrema_values[1] - lambda_vec[2:llv_1])))

  # Get new pdf area
  area_pdf_diff <- abs(trapz(X, exp(log_like)) - 1)

  diffeo_list = list(number_extrema = found_extrema$num_extrema,
                     extrema_heights_diff = extrema_heights_diff,
                     area_pdf_diff = area_pdf_diff)

  return(diffeo_list)

}

diffeo_list_1 <- test_diffeo_log_like(X, Beta_1, b_vec_1, lambda_vec_1)
diffeo_list_2 <- test_diffeo_log_like(X, Beta_2, b_vec_2, lambda_vec_2)
diffeo_list_3 <- test_diffeo_log_like(X, Beta_3, b_vec_3, lambda_vec_3)

num_extrema_1 <- length(b_vec_1) - 2
num_extrema_2 <- length(b_vec_2) - 2
num_extrema_3 <- length(b_vec_3) - 2

test_that('Test if number of extrema remains the same',
          {expect_equal(diffeo_list_1$number_extrema, num_extrema_1)
            expect_equal(diffeo_list_2$number_extrema, num_extrema_2)
            expect_equal(diffeo_list_3$number_extrema, num_extrema_3)})

test_that('Test if height ratios remains the same',
          {expect_lte(diffeo_list_1$extrema_heights_diff, 1e-4)
            expect_lte(diffeo_list_2$extrema_heights_diff, 1e-4)
            expect_lte(diffeo_list_3$extrema_heights_diff, 1e-4)})

test_that('Test if pdf area is approximately 1 (error higher because of trapezoid integration)',
          {expect_lte(diffeo_list_1$area_pdf_diff, 0.05)
            expect_lte(diffeo_list_2$area_pdf_diff, 0.05)
            expect_lte(diffeo_list_3$area_pdf_diff, 0.05)})



# Initialize Data
set.seed(42)
X1 <- rnorm(75)
X2 <- rgamma(50, 2, 1)
X3 <- c(runif(75), rnorm(25, 0.75, 0.001))

# Initialize set beta_starts
beta_starts_1 <- matrix(runif(25 * 3), ncol = 3)
beta_starts_2 <- matrix(runif(25 * 2), ncol = 2)
beta_starts_3 <- matrix(runif(25 * 4), ncol = 4)

# Find MLE and see if close to truth
mle_1 <- diffeo_mle_estimate(X1, 3, beta_starts = beta_starts_1,
                             plot_results = FALSE)
mle_2 <- diffeo_mle_estimate(X2, 2, beta_starts = beta_starts_2,
                             plot_results = FALSE)
mle_3 <- diffeo_mle_estimate(X3, 4, beta_starts = beta_starts_3,
                             plot_results = FALSE)

test_that('MLE produces proper mode estimates',
          {expect_lte(abs(mle_1$mode_estimate), 0.3)
            expect_lte(abs(mle_2$mode_estimate - 1), 0.3)
            expect_lte(abs(mle_3$mode_estimate - 0.75), 0.1)})

# Find Bayes estimates and see if close to truth
bayes_1 <- diffeo_bayes_estimate(X1, 3, beta_starts = beta_starts_1,
                             plot_results = FALSE)
bayes_2 <- diffeo_bayes_estimate(X2, 2, beta_starts = beta_starts_2,
                             plot_results = FALSE)
bayes_3 <- diffeo_bayes_estimate(X3, 4, beta_starts = beta_starts_3,
                             plot_results = FALSE)

# MAP test
test_that('MAP produces proper mode estimates',
          {expect_lte(abs(bayes_1$bayes_map_mode), 0.3)
            expect_lte(abs(bayes_2$bayes_map_mode - 1), 0.3)
            expect_lte(abs(bayes_3$bayes_map_mode - 0.75), 0.1)})

# Test if true mode is in 95% credible interval
quantile_1 <- quantile(bayes_1$sampled_modes, probs = c(0.025, 0.975))
quantile_2 <- quantile(bayes_2$sampled_modes, probs = c(0.025, 0.975))
quantile_3 <- quantile(bayes_3$sampled_modes, probs = c(0.025, 0.975))

test_that('MAP produces proper mode estimates',
          {expect_lte(as.numeric(quantile_1[1]), 0)
            expect_gte(as.numeric(quantile_1[2]), 0)
            expect_lte(as.numeric(quantile_2[1]), 1)
            expect_gte(as.numeric(quantile_2[2]), 1)
            expect_lte(as.numeric(quantile_3[1]), 0.75)
            expect_gte(as.numeric(quantile_3[2]), 0.75)})

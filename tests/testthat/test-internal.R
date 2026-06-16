# The Dirichlet density is inlined (ddirichlet_local) to drop the gtools dependency.
# Verify it equals the closed-form Dirichlet pdf computed independently.

test_that("ddirichlet_local matches the closed-form Dirichlet density", {
  ref <- function(x, alpha) {
    exp(lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(x)))
  }
  set.seed(42)
  for (i in seq_len(500)) {
    x <- runif(2); x <- x / sum(x)
    a <- runif(2, 0.1, 50)
    expect_equal(rabhit:::ddirichlet_local(x, a), ref(x, a))
  }
})

test_that("ddirichlet_local returns 0 for invalid input (matching gtools semantics)", {
  expect_identical(rabhit:::ddirichlet_local(c(0.6, 0.6), c(1, 1)), 0) # does not sum to 1
  expect_identical(rabhit:::ddirichlet_local(c(-0.2, 1.2), c(1, 1)), 0) # outside [0, 1]
})

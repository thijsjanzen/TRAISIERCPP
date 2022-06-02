context("sample_event")

test_that("sample_event use", {
  v <- 1:10
  num_repl <- 1e3
  a1 <- sample(v, size = num_repl, replace = TRUE)
  ax <- table(a1)
  testthat::expect_equal(ax[[1]], num_repl*0.1, tolerance = 100)

  a2 <- c()
  for (r in 1:1000) {
    a2[r] <- TRAISIERCPP::sample_event(rep(0.1, 10))
  }
  ax <- table(a2)
  testthat::expect_equal(ax[[1]], num_repl*0.1, tolerance = 100)

  a2 <- c()
  probs <- rep(0.05, 10)
  probs[1] <- 10
  for (r in 1:1000) {
    a2[r] <- TRAISIERCPP::sample_event(probs)
  }
  ax <- table(a2)
  ax <- ax / sum(ax)
  for (i in 1:length(probs)) {
    testthat::expect_equal(ax[[i]], probs[i] / sum(probs), tolerance = 0.05)
  }
})

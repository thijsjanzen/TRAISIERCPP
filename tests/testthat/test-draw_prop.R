context("test_draw_prop")

test_that("draw_prop", {

  probs <- c(0.1, 0.8, 0.05, 0.05)
  found <- rep(NA, 1000)
  for (r in 1:1000) {
    found[r] <- TRAISIERCPP::test_draw_prop(probs)
  }
  a1 <- table(found)
  ax <- a1 / sum(a1)
  for (i in 1:length(probs)) {
    testthat::expect_equal(ax[[i]], probs[i], tolerance = 0.03)
  }
})

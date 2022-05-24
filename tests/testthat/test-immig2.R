context("immigration2")

test_that("immigration2", {
  # first, we test sample_spec

  get_R_colonist <- function(mainland_spec_, M2_) {
    mainland1 = length(mainland_spec_)
    mainland2 = M2_
    mainland_total = mainland1 + mainland2
    return(DDD::sample2((mainland1 + 1):mainland_total, 1))
  }


  mainland_n <- 10

  mainland_spec <- seq(1, mainland_n, 1)
  M2 <- 10
  found_R <- c()
  found_C <- c()
  for (r in 1:10000) {
    found_R[r] <- get_R_colonist(mainland_spec, M2)
    found_C[r] <- TRAISIERCPP::test_sample_spec(mainland_spec, M2)
  }

  a1 <- table(found_R)
  a1 <- a1 / sum(a1)
  a2 <- table(found_C)
  a2 <- a2 / sum(a2)
  testthat::expect_equal(names(a1), names(a2))
  for (i in 1:length(a1)) {
    testthat::expect_equal(a1[[i]], a2[[i]], tolerance = 0.02)
  }


  # and now we do actual immigration

  timeval = 1.42
  mainland_spec <- c(1)
  island_spec <- rbind(rep("1", 8),
                       rep("1", 8))
  island_spec[, 4] <- "I"

  res <- TRAISIERCPP::test_immigration2(island_spec,
                                       mainland_spec,
                                       timeval,
                                       2)

  testthat::expect_true(length(res[, 1]) == 3)
  testthat::expect_true(as.numeric(res[3, 3]) == timeval)
  testthat::expect_true(res[3, 4] == "I")
  testthat::expect_equal(as.numeric(res[3, 8]), 2)
})

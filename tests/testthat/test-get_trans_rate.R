context("get_trans_rate")

test_that("get_trans_rate", {
  trans_rate1 <- 0.01
  trans_rate2 <- 0.02

  trait_pars <- c()
  trait_pars$trans_rate <- trans_rate1
  trait_pars$trans_rate2 <- trans_rate2

  island_spec <- rbind(rep("1", 8),
                       rep("1", 8),
                       rep("2", 8),
                       rep("2", 8),
                       rep("2", 8))
  island_spec[, 4] <- "I"

  # all are immig, so no need to check for column 4
  num_spec_trait1 <- length(which(island_spec[, 8] == "1"))
  num_spec_trait2 <- length(which(island_spec[, 8] == "2"))


  reference_value <- TRAISIERCPP::DAISIE_get_trans_rate(trait_pars, island_spec)

  rcpp_value <- TRAISIERCPP::test_get_trans_rate(trans_rate1,
                                                 trans_rate2,
                                                 num_spec_trait1,
                                                 num_spec_trait2)

  testthat::expect_equal(reference_value$trans_rate1,
                         rcpp_value[[1]])
  testthat::expect_equal(reference_value$trans_rate2,
                         rcpp_value[[2]])

})

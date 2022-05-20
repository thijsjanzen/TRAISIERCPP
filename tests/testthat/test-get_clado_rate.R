context("get_clado_rate")

test_that("get_clado_rate", {
  lac <- 0.01
  num_spec <- 10
  A <- 1
  clado_rate2 <- 0.02

  trait_pars <- c()
  trait_pars$clado_rate2 <- clado_rate2

  island_spec <- rbind(rep("1", 8),
                       rep("1", 8),
                       rep("2", 8),
                       rep("2", 8),
                       rep("2", 8))
  island_spec[, 4] <- "I"

  # all are immig, so no need to check for column 4
  num_spec_trait1 <- length(which(island_spec[, 8] == "1"))
  num_spec_trait2 <- length(which(island_spec[, 8] == "2"))


  for (K in c(10, 100, Inf)) {

    reference_value <- TRAISIERCPP::DAISIE_get_clado_rate(lac,
                                                          hyper_pars = TRAISIERCPP::create_hyper_pars(d = 0, x = 0),
                                                          num_spec,
                                                          K,
                                                          A,
                                                          trait_pars,
                                                          island_spec)

    rcpp_value <- TRAISIERCPP::test_get_clado_rate(lac,
                                                   num_spec,
                                                   K,
                                                   A,
                                                   clado_rate2,
                                                   num_spec_trait1,
                                                   num_spec_trait2)

    testthat::expect_equal(reference_value$clado_rate1,
                           rcpp_value[[1]])
    testthat::expect_equal(reference_value$clado_rate2,
                           rcpp_value[[2]])
  }
})

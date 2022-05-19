context("get_immig_rate")

test_that("get_immig_rate", {

  M2 <- 1
  immig_rate2 <- 0.2
  trait_pars <- c()
  trait_pars$M2 <- M2
  trait_pars$immig_rate2 <- immig_rate2

  gam <- 0.1
  A <- 1
  num_spec <- 10

  mainland_n <- 1

  for (K in c(10, 100, Inf)) {

    reference_value <- TRAISIERCPP::DAISIE_get_immig_rate(gam,
                                                          A,
                                                          num_spec,
                                                          K,
                                                          mainland_n,
                                                          trait_pars)

    rcpp_value <- TRAISIERCPP::test_get_immig_rate(gam,
                                                   A,
                                                   num_spec,
                                                   K,
                                                   mainland_n,
                                                   M2,
                                                   immig_rate2)

    testthat::expect_equal(reference_value$immig_rate1,
                           rcpp_value[[1]])
    testthat::expect_equal(reference_value$immig_rate2,
                           rcpp_value[[2]])
  }
})

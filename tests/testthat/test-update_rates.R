context("update_rates")

test_that("update_rates", {
  trait_pars <- TRAISIERCPP::create_trait_pars(trans_rate = 0.01,
                                               immig_rate2 = 0.02,
                                               ext_rate2 = 0.02,
                                               ana_rate2 = 0.02,
                                               clado_rate2 = 0.02,
                                               trans_rate2 = 0.02,
                                               M2 = 1)

  timeval <- 0.01
  total_time <- 2

  gam = immig_rate1 = 0.01
  laa = ana_rate1 = 0.01
  lac = clado_rate1 <- 0.01
  mu = ext_rate1 <- 0.01
  extcutoff = 1000
  K <- Inf
  num_spec <- 10
  num_immigrants <- 10
  mainland_n <- 1


  island_spec <- rbind(rep("1", 8),
                       rep("1", 8),
                       rep("2", 8),
                       rep("2", 8),
                       rep("2", 8))
  island_spec[, 4] <- "I"



  reference_value <- TRAISIERCPP::DAISIE_update_rates_trait(timeval = timeval,
                                                            total_time = total_time,
                                                            gam = gam,
                                                            laa = laa,
                                                            lac = lac,
                                                            mu = mu,
                                                            hyper_pars = TRAISIERCPP::create_hyper_pars(d = 0, x = 0),
                                                            extcutoff = extcutoff,
                                                            K = K,
                                                            num_spec = num_spec,
                                                            num_immigrants = num_immigrants,
                                                            mainland_n = mainland_n,
                                                            trait_pars = trait_pars,
                                                            island_spec = island_spec)


  # all are immig, so no need to check for column 4

  island_spec_R <- rbind(rep(1, 8),
                       rep(1, 8),
                       rep(2, 8),
                       rep(2, 8),
                       rep(2, 8))

  rcpp_value <- TRAISIERCPP::test_update_rates(timeval,
                                               total_time,
                                               gam,
                                               laa,
                                               lac,
                                               mu,
                                               K,
                                               num_spec,
                                               num_immigrants,
                                               mainland_n,
                                               island_spec_R,
                                               trait_pars)


  testthat::expect_equal(reference_value$immig_rate,
                         rcpp_value[[1]])
  testthat::expect_equal(reference_value$ext_rate,
                         rcpp_value[[2]])
  testthat::expect_equal(reference_value$ana_rate,
                         rcpp_value[[3]])
  testthat::expect_equal(reference_value$clado_rate,
                         rcpp_value[[4]])
  testthat::expect_equal(reference_value$trans_rate,
                         rcpp_value[[5]])

  testthat::expect_equal(reference_value$immig_rate2,
                         rcpp_value[[6]])
  testthat::expect_equal(reference_value$ext_rate2,
                         rcpp_value[[7]])
  testthat::expect_equal(reference_value$ana_rate2,
                         rcpp_value[[8]])
  testthat::expect_equal(reference_value$clado_rate2,
                         rcpp_value[[9]])
  testthat::expect_equal(reference_value$trans_rate2,
                         rcpp_value[[10]])

  testthat::expect_equal(reference_value$M2,
                         rcpp_value[[11]])

})


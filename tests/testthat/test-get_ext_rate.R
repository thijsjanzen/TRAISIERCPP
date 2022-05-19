context("get_immig_rate")

test_that("get_immig_rate", {
  mu <- 0.01
  num_spec <- 5
  A <- 1
  ext_rate2 <- 0.02
  num_spec_trait1 <- 5
  num_spec_trait2 <- 10

  trait_pars <- c()
  trait_pars$ext_rate2 <- ext_rate2

  island_spec <- rbind(rep("1", 8),
                       rep("1", 8),
                       rep("2", 8),
  rep("2", 8),
  rep("2", 8))

num_spec_trait1 <- length(which(island_spec[, 8] == "1"))
num_spec_trait2 <- length(which(island_spec[, 8] == "2"))


reference_value <- TRAISIERCPP::DAISIE_get_ext_rate(mu,
                        hyper_pars = DAISIE::create_hyper_pars(d = 0, x = 0),
                        num_spec =  num_spec,
                                                    A =          A,
                                                    island_spec = island_spec,
                                                    trait_pars = trait_pars)

rcpp_value <- TRAISIERCPP::test_get_ext_rate(mu,
                                             num_spec,
                                             A,
                                             ext_rate2,
                                             num_spec_trait1,
                                             num_spec_trait2)

testthat::expect_equal(reference_value$ext_rate1,
                       rcpp_value[[1]])
testthat::expect_equal(reference_value$ext_rate2,
                       rcpp_value[[2]])
})

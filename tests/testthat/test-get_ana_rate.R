context("get_ana_rate")

test_that("get_ana_rate", {
  laa <- 0.01
  num_immigrants <- 10
  ana_rate2 <- 0.02

  trait_pars <- c()
  trait_pars$ana_rate2 <- ana_rate2

  island_spec <- rbind(rep("1", 8),
                       rep("1", 8),
                       rep("2", 8),
  rep("2", 8),
  rep("2", 8))
  island_spec[, 4] <- "I"

  # all are immig, so no need to check for column 4
num_immig_trait1 <- length(which(island_spec[, 8] == "1"))
num_immig_trait2 <- length(which(island_spec[, 8] == "2"))


reference_value <- TRAISIERCPP::DAISIE_get_ana_rate(laa,
                                                    num_immigrants,
                                                    island_spec,
                                                    trait_pars)

rcpp_value <- TRAISIERCPP::test_get_ana_rate(laa,
                                             num_immigrants,
                                             ana_rate2,
                                             num_immig_trait1,
                                             num_immig_trait2)

testthat::expect_equal(reference_value$ana_rate1,
                       rcpp_value[[1]])
testthat::expect_equal(reference_value$ana_rate2,
                       rcpp_value[[2]])
})

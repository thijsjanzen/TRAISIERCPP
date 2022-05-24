context("cladogenesis")

test_that("cladogenesis", {
  set.seed(1)
  island_spec_R <- TRAISIERCPP::create_island_spec(time = 3,
                                                   mainland_n = 3,
                                                   K = Inf,
                                                   immig_rate = 0.1,
                                                   ext_rate = 0.0,
                                                   ana_rate = 0.0,
                                                   clado_rate = 0.2,
                                                   immig_rate2 = 0.0,
                                                   ext_rate2 = 0.0,
                                                   ana_rate2 = 0.0,
                                                   clado_rate2 = 0.0,
                                                   trans_rate = 0.0,
                                                   trans_rate2 = 0.0,
                                                   M2 = 1)

  island_spec_R2 <- island_spec_R
  island_spec_R2[is.na(island_spec_R2[, 5]), 5] <- "D"

  island_spec_R2[2, 8] <- "2"

  for (focal_trait in c(1, 2)) {
    island_spec_R3 <- island_spec_R2
    island_spec_R3 <- test_transition(island_spec_R3, focal_trait)
    num_trait <- as.numeric(island_spec_R3[, 8])
    testthat::expect_equal(sum(num_trait == focal_trait), 0)
  }
})

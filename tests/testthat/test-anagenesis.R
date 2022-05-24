context("anagenesis")

test_that("anagenesis", {
  set.seed(1)
  island_spec_R <- TRAISIERCPP::create_island_spec(time = 3,
                                                   mainland_n = 3,
                                                   K = Inf,
                                                   immig_rate = 0.8,
                                                   ext_rate = 0.0,
                                                   ana_rate = 0.0,
                                                   clado_rate = 0.0,
                                                   immig_rate2 = 0.2,
                                                   ext_rate2 = 0.0,
                                                   ana_rate2 = 0.0,
                                                   clado_rate2 = 0.0,
                                                   trans_rate = 0.0,
                                                   trans_rate2 = 0.0,
                                                   M2 = 1)


  for (focal_trait in 1:2) {
    island_spec_R2 <- island_spec_R
    island_spec_R2[is.na(island_spec_R2[, 5]), 5] <- "D"
    if (focal_trait == 1) {
       island_spec_R2[,  8] <- "2"
       island_spec_R2[1, 8] <- "1"
    } else {
      island_spec_R2[, 8] <- "1"
      island_spec_R2[1, 8] <- "2"
    }
    max_id <- 1 + max(as.numeric(island_spec_R2[, 1]))
    res <- TRAISIERCPP::test_anagenesis(island_spec_R2, max_id, focal_trait)

    testthat::expect_equal(res[1, 4], "A")
    testthat::expect_equal(res[1, 7], "Immig_parent")
    testthat::expect_equal(as.numeric(res[1, 8]), focal_trait)
  }
})


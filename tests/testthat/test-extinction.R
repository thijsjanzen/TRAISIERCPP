context("extinction")

test_that("extinction full", {
  set.seed(1)
  island_spec_R <- TRAISIERCPP::create_island_spec(time = 5,
                                                   mainland_n = 10,
                                                   K = Inf,
                                                   immig_rate = 0.2,
                                                   ext_rate = 0.0,
                                                   ana_rate = 0.2,
                                                   clado_rate = 0.2,
                                                   immig_rate2 = 0.1,
                                                   ext_rate2 = 0.0,
                                                   ana_rate2 = 0.1,
                                                   clado_rate2 = 0.1,
                                                   trans_rate = 0.1,
                                                   trans_rate2 = 0.1,
                                                   M2 = 1)

  island_spec_R2 <- island_spec_R
  island_spec_R2[is.na(island_spec_R2[, 5]), 5] <- "D"

  for (index in 1:length(island_spec_R[, 1])) {
    table1 <- TRAISIERCPP::DAISIE_test_execute_extinction(island_spec_R,
                                                          index)
    table2 <- TRAISIERCPP::test_execute_extinction(island_spec_R2,
                                                   index)
    testthat::expect_equal(length(table1[, 1]),
                           length(table2[, 1]))
    for (i in 1:length(table1[, 1])) {
      testthat::expect_equal(as.numeric(table1[i, 1]),
                             as.numeric(table2[i, 1]))

      testthat::expect_equal(as.numeric(table1[i, 3]),
                             as.numeric(table2[i, 3]), tolerance = 0.1)

      testthat::expect_equal(table1[i, 4],
                             table2[i, 4])

      testthat::expect_equal(as.numeric(table1[i, 8]),
                             as.numeric(table2[i, 8]))
    }
  }
})



test_that("extinction simple" , {

  island_spec_R <- rbind(rep("1", 8),
                         rep("1", 8),
                         rep("1", 8),
                         rep("2", 8),
                         rep("2", 8),
                         rep("2", 8))

  num_trait_1 <- length(island_spec_R[, 8] == "1")
  num_trait_2 <- length(island_spec_R[, 8] == "2")
  island_spec_R[, 4] <- "I"

  for (focal_trait in c(1, 2)) {

    res <- TRAISIERCPP::test_extinction(island_spec_R, focal_trait)

    if (focal_trait == 1) {
      num_trait_1_res <- length(res[, 8] == "1")
      testthat::expect_true(num_trait_1 > num_trait_1_res)
    }
    if (focal_trait == 2) {
      num_trait_2_res <- length(res[, 8] == "2")
      testthat::expect_true(num_trait_2 > num_trait_2_res)
    }
  }
})

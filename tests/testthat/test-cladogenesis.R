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

  # first, we test for type C

    island_spec_R2 <- island_spec_R
    island_spec_R2[is.na(island_spec_R2[, 5]), 5] <- "D"

    max_id <- 1 + max(as.numeric(island_spec_R2[, 1]))
    timeval <- 3
    focal_trait <- 1
    island_spec_R2[2, 8] <- "2"

    res <- TRAISIERCPP::test_cladogenesis(island_spec_R2, timeval,
                                          max_id, focal_trait)

    testthat::expect_equal(length(res[, 1]), length(island_spec_R2[, 1]) + 1)
    testthat::expect_equal(as.numeric(res[1, 2]), max_id + 1)
    testthat::expect_equal(as.numeric(res[3, 1]), max_id + 2)
    testthat::expect_equal(as.numeric(res[1, 3]), as.numeric(res[3,3]))
    testthat::expect_equal(res[3, 4], "C")
    tail_spec <- res[1, 5]
    tail_spec <- substr(tail_spec, start = 2, stop = 2)
    testthat::expect_equal(tail_spec, "A")

    tail_spec <- res[3, 5]
    tail_spec <- substr(tail_spec, start = 2, stop = 2)
    testthat::expect_equal(tail_spec, "B")


    # now, we test for type A

    island_spec_R2 <- island_spec_R
    island_spec_R2[is.na(island_spec_R2[, 5]), 5] <- "D"
    island_spec_R2[1, 4] <- "A"

    max_id <- 1 + max(as.numeric(island_spec_R2[, 1]))
    timeval <- 3
    focal_trait <- 1
    island_spec_R2[2, 8] <- "2"

    res <- TRAISIERCPP::test_cladogenesis(island_spec_R2, timeval,
                                          max_id, focal_trait)

    testthat::expect_equal(res[1, 4], "C")
    testthat::expect_equal(res[3, 4], "C")

    testthat::expect_equal(as.numeric(res[1, 3]), as.numeric(res[1, 6]))

    testthat::expect_equal(length(res[, 1]), length(island_spec_R2[, 1]) + 1)

    testthat::expect_equal(as.numeric(res[1, 2]), max_id + 1)
    testthat::expect_equal(as.numeric(res[3, 1]), max_id + 2)

    testthat::expect_equal(as.numeric(res[1, 3]), as.numeric(res[3, 3]))

    testthat::expect_equal(as.numeric(res[3, 8]), focal_trait)

    tail_spec <- res[1, 5]
    testthat::expect_equal(tail_spec, "A")

    tail_spec <- res[3, 5]
    testthat::expect_equal(tail_spec, "B")
})


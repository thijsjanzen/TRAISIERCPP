context("immigration")

test_that("immigration", {
   timeval = 1.42
   mainland_spec <- c(1)
   island_spec <- rep("1", 8)
   island_spec[4] <- "I"

   res <- TRAISIERCPP::test_immigration(timeval,
                                        mainland_spec,
                                        island_spec)
   testthat::expect_true(length(res[, 1]) == 2)
   testthat::expect_true(as.numeric(res[2, 3]) == timeval)
})

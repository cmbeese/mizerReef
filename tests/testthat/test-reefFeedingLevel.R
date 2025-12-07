test_that("reefFeedingLevel runs without error", {
    params <- newReefParams()
    expect_error(reefFeedingLevel(params, n = NULL, n_pp = NULL, n_other = NULL, t = 0, encounter = NULL), NA)
})

test_that("reefPredMort runs without error", {
    params <- newReefParams()
    expect_error(reefPredMort(params, n = NULL, n_pp = NULL, n_other = NULL, t = 0, pred_rate = NULL), NA)
})

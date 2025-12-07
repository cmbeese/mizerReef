test_that("reefMort runs without error", {
    params <- newReefParams()
    expect_error(reefMort(params, n = NULL, n_pp = NULL, n_other = NULL, t = 0, f_mort = NULL, pred_mort = NULL), NA)
})

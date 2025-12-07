test_that("reefRates runs without error", {
    params <- newReefParams()
    expect_error(reefRates(params, n = NULL, n_pp = NULL, n_other = NULL, t = 0, effort = NULL, rates_fns = NULL), NA)
})

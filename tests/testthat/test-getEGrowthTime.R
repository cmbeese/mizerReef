test_that("getEGrowthTime runs without error", {
    params <- newReefParams()
    expect_error(getEGrowthTime(params, n = NULL, n_pp = NULL, n_other = NULL), NA)
})

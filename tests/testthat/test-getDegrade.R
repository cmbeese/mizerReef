test_that("getDegrade runs without error", {
    params <- newReefParams()
    expect_error(getDegrade(params, n = NULL, n_pp = NULL, n_other = NULL), NA)
})

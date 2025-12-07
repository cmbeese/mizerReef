test_that("getVulnerable runs without error", {
    params <- newReefParams()
    expect_error(getVulnerable(params, n = NULL, n_pp = NULL, n_other = NULL), NA)
})

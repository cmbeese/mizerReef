test_that("reefVulnerable runs without error", {
    params <- newReefParams()
    expect_error(reefVulnerable(params, n = NULL, n_pp = NULL, n_other = NULL, t = 0), NA)
})

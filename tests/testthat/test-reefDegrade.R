test_that("reefDegrade runs without error", {
    params <- newReefParams()
    expect_s4_class(reefDegrade(params, n = NULL, n_pp = NULL, n_other = NULL, t = 1), "MizerReefParams")
})

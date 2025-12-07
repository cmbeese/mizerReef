test_that("reefEncounter runs without error", {
    params <- newReefParams()
    expect_error(reefEncounter(params, n = NULL, n_pp = NULL, n_other = NULL, t = 0), NA)
})

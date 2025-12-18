test_that("getVulnerable runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(getVulnerable(params, 
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other), NA)
})

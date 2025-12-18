test_that("reefPredMort runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    r <- getRates(params)
    expect_error(reefPredMort(params, 
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other,
        t = 0,
        pred_rate = r$pred_rate,
        vulnerable = r$vulnerable), NA)
})
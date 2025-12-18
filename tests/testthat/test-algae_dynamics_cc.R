test_that("algae_dynamics_cc runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(algae_dynamics_cc(params,
        n = params@initial_n,
        n_other = params@initial_n_other,
        rates = getRates(params),
        dt = 1), NA)
})

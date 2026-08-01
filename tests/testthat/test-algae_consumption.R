test_that("algae_consumption matches an independently computed sum", {
    # Independent computation via per-species summation (not the source
    # function's matrix-multiply) of rho * n * dw.
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    rho <- params@other_params$algae$rho

    expected <- 0
    for (sp in seq_len(nrow(n))) {
        expected <- expected + sum(rho[sp, ] * n[sp, ] * params@dw)
    }

    expect_equal(algae_consumption(params, n = n, rates = getRates(params)),
                 expected)
})

test_that("algae_consumption ignores feeding level, unlike detritus_consumption", {
    # Documented design decision (see the "Algae consumption" roxygen
    # section): algal grazing pressure is continuous and does not depend on
    # individual consumer satiation. Confirm the rates argument's
    # feeding_level truly has no effect on the result.
    data(caribbean_3_model)
    params <- caribbean_3_model
    real_rates <- getRates(params)
    fake_rates <- real_rates
    fake_rates$feeding_level[] <- 0.9999

    expect_equal(
        algae_consumption(params, rates = real_rates),
        algae_consumption(params, rates = fake_rates)
    )
})

test_that("algae_consumption is zero when no species interacts with algae", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@other_params$algae$rho[] <- 0
    expect_equal(algae_consumption(params), 0)
})

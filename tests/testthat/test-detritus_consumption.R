test_that("detritus_consumption matches an independently computed sum with feeding level applied", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    rates <- getRates(params)
    rho <- params@other_params$detritus$rho

    expected <- 0
    for (sp in seq_len(nrow(n))) {
        expected <- expected + sum(rho[sp, ] * n[sp, ] *
            (1 - rates$feeding_level[sp, ]) * params@dw)
    }

    expect_equal(detritus_consumption(params, n = n, rates = rates), expected)
})

test_that("detritus_consumption drops to zero when feeding level is fully satiated", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    rates <- getRates(params)
    rates$feeding_level[] <- 1

    expect_equal(detritus_consumption(params, rates = rates), 0)
})

test_that("detritus_consumption equals the undiscounted rho*n*dw sum when feeding level is zero", {
    # Contrast with algae_consumption(), which is *always* undiscounted --
    # this confirms detritus_consumption's feeding-level factor is a real,
    # switchable effect rather than a no-op.
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    rates <- getRates(params)
    rates$feeding_level[] <- 0

    expected <- sum((params@other_params$detritus$rho * n) %*% params@dw)
    expect_equal(detritus_consumption(params, n = n, rates = rates), expected)
})

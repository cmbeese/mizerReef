test_that("getDetritusProduction's external component is exactly other_params$detritus$external", {
    # Regression guard for the fix that consolidated detritus params onto
    # the single canonical other_params$detritus location: this is exactly
    # the field that used to go stale.
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@other_params$detritus$external <- 12345
    expect_equal(unname(getDetritusProduction(params)["external"]), 12345)
})

test_that("getDetritusProduction's feces component matches an independent per-species computation", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    rates <- getRates(params)
    alpha <- params@species_params$alpha

    expected_feces <- 0
    for (sp in seq_len(nrow(n))) {
        consumed <- (1 - rates$feeding_level[sp, ]) * rates$encounter[sp, ] *
            n[sp, ] * params@dw
        expected_feces <- expected_feces + sum((1 - alpha[sp]) * consumed)
    }

    expect_equal(unname(getDetritusProduction(params, n, rates)["feces"]),
                 expected_feces)
})

test_that("getDetritusProduction's decomp component correctly weights senescence vs residual mortality", {
    # Uses the already-independently-tested reefSenMort()/params@mu_b
    # building blocks directly (not re-derived) to check that
    # getDetritusProduction() combines them with the documented
    # sen_decomp/ext_decomp weights and doesn't, e.g., have them swapped.
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n

    ex_mort <- sum((params@mu_b * n) %*% (params@w * params@dw))
    sen_mort <- sum((getSenMort(params) * n) %*% (params@w * params@dw))
    sen_decomp <- params@other_params$detritus$sen_decomp
    ext_decomp <- params@other_params$detritus$ext_decomp
    expected_decomp <- (ext_decomp * ex_mort) + (sen_decomp * sen_mort)

    expect_equal(unname(getDetritusProduction(params)["decomp"]), expected_decomp)

    # Confirm the weights are not swapped: using distinct sen_decomp/
    # ext_decomp values should change the result in the direction implied
    # by which mortality source dominates.
    skip_if(sen_mort == ex_mort, "sen_mort and ex_mort coincide, can't distinguish weighting")
    params@other_params$detritus$sen_decomp <- 1
    params@other_params$detritus$ext_decomp <- 0
    expect_equal(unname(getDetritusProduction(params)["decomp"]), sen_mort)
})

test_that("getDetritusProduction returns named feces/decomp/external entries", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_named(getDetritusProduction(params), c("feces", "decomp", "external"))
})

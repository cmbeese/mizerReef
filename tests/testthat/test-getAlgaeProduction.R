test_that("getAlgaeProduction runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(getAlgaeProduction(params), NA)
})

test_that("getAlgaeProduction applies the growth boost multiplier to the baseline rate", {
    data(caribbean_3_model)
    data(rubble_scale)
    growth_boost <- c(1.11, 1.22, 1.33, 1.11)
    params <- setDegradation(caribbean_3_model,
        trajectory = "rubble", deg_scale = rubble_scale,
        bleach_time = 2, degrade = TRUE, algae_boost = TRUE,
        algae_growth_boost = growth_boost, algae_capacity_boost = c(2)
    )
    baseline <- params@other_params$algae$growth

    expect_equal(getAlgaeProduction(params, t = 0), baseline)
    expect_equal(getAlgaeProduction(params, t = 2), baseline * growth_boost[1])
    expect_equal(
        getAlgaeProduction(params, t = 3),
        baseline * growth_boost[1] * growth_boost[2]
    )
})

test_that("getAlgaeProduction defaults to t = 0 (no boost) when t is not given", {
    data(caribbean_3_model)
    data(rubble_scale)
    params <- setDegradation(caribbean_3_model,
        trajectory = "rubble", deg_scale = rubble_scale,
        bleach_time = 2, degrade = TRUE, algae_boost = TRUE,
        algae_growth_boost = c(1.5)
    )
    expect_equal(
        getAlgaeProduction(params),
        params@other_params$algae$growth
    )
})

test_that("algae_boost = FALSE leaves algae production at the baseline rate regardless of t", {
    data(caribbean_3_model)
    data(rubble_scale)
    params <- setDegradation(caribbean_3_model,
        trajectory = "rubble", deg_scale = rubble_scale,
        bleach_time = 2, degrade = TRUE, algae_boost = FALSE,
        algae_growth_boost = c(5), algae_capacity_boost = c(5)
    )
    baseline <- params@other_params$algae$growth
    expect_equal(getAlgaeProduction(params, t = 2), baseline)
    expect_equal(getAlgaeProduction(params, t = 10), baseline)
})

test_that("getAlgaeBoost returns 1 (no boost) when algae_boost is not set", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_identical(
        getAlgaeBoost(params, t = 5, boost_vector = c(1.5, 2)),
        1
    )
})

test_that("getAlgaeBoost returns 1 before the bleaching year", {
    data(caribbean_3_model)
    data(rubble_scale)
    params <- setDegradation(caribbean_3_model,
        trajectory = "rubble", deg_scale = rubble_scale,
        bleach_time = 3, degrade = TRUE, algae_boost = TRUE,
        algae_growth_boost = c(1.5, 2), algae_capacity_boost = c(2)
    )
    boost_vec <- params@other_params$algae_params$algae_growth_boost
    expect_identical(getAlgaeBoost(params, t = 0, boost_vec), 1)
    expect_identical(getAlgaeBoost(params, t = 2, boost_vec), 1)
})

test_that("getAlgaeBoost compounds across post-bleaching years and then plateaus", {
    data(caribbean_3_model)
    data(rubble_scale)
    growth_boost <- c(1.11, 1.22, 1.33, 1.11)
    params <- setDegradation(caribbean_3_model,
        trajectory = "rubble", deg_scale = rubble_scale,
        bleach_time = 2, degrade = TRUE, algae_boost = TRUE,
        algae_growth_boost = growth_boost, algae_capacity_boost = c(2)
    )

    # At the bleach year, only the first element applies
    expect_equal(getAlgaeBoost(params, t = 2, growth_boost), growth_boost[1])
    # Each subsequent year multiplies in the next element (cumulative product)
    expect_equal(
        getAlgaeBoost(params, t = 3, growth_boost),
        growth_boost[1] * growth_boost[2]
    )
    expect_equal(
        getAlgaeBoost(params, t = 4, growth_boost),
        growth_boost[1] * growth_boost[2] * growth_boost[3]
    )
    full_product <- prod(growth_boost)
    expect_equal(getAlgaeBoost(params, t = 5, growth_boost), full_product)
    # Once t runs past the length of the boost vector, the cumulative
    # product plateaus rather than continuing to compound
    expect_equal(getAlgaeBoost(params, t = 6, growth_boost), full_product)
    expect_equal(getAlgaeBoost(params, t = 20, growth_boost), full_product)
})

test_that("getAlgaeBoost returns 1 for an empty or NULL boost vector", {
    data(caribbean_3_model)
    data(rubble_scale)
    params <- setDegradation(caribbean_3_model,
        trajectory = "rubble", deg_scale = rubble_scale,
        bleach_time = 2, degrade = TRUE, algae_boost = TRUE
    )
    expect_identical(getAlgaeBoost(params, t = 5, boost_vector = NULL), 1)
    expect_identical(getAlgaeBoost(params, t = 5, boost_vector = numeric(0)), 1)
})

test_that("a full simulation with degrade and algae_boost both TRUE runs without error", {
    data(caribbean_3_model)
    data(rubble_scale)
    params <- setDegradation(caribbean_3_model,
        trajectory = "rubble", deg_scale = rubble_scale,
        bleach_time = 2, degrade = TRUE, algae_boost = TRUE,
        algae_growth_boost = c(1.11, 1.22, 1.33, 1.11),
        algae_capacity_boost = c(2.0, 1.5)
    )
    expect_error(project(params, t_max = 6, progress_bar = FALSE), NA)
})

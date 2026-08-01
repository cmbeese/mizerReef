test_that("calibrateReefNumber returns params unchanged when number_observed is absent", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_false("number_observed" %in% names(params@species_params))
    expect_equal(calibrateReefNumber(params), params)
})

test_that("calibrateReefNumber rescales by exactly observed_total / model_total", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    no_sp <- nrow(params@species_params)
    params@species_params$number_observed <- c(1000, NA, 500)[seq_len(no_sp)]
    params@species_params$number_cutoff <- rep(0, no_sp)

    observed <- params@species_params$number_observed
    observed_total <- sum(observed, na.rm = TRUE)
    model_total <- 0
    for (sp_idx in which(!is.na(observed))) {
        model_total <- model_total + sum(params@initial_n[sp_idx, ] * params@dw)
    }
    expected_factor <- observed_total / model_total

    expect_equal(calibrateReefNumber(params), scaleReefModel(params, factor = expected_factor))
})

test_that("calibrateReefNumber returns params unchanged when number_observed is all NA", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@species_params$number_observed <- NA
    expect_equal(calibrateReefNumber(params), params)
})

test_that("calibrateReefNumber makes total model number match total observed", {
    # End-to-end postcondition check, mirroring calibrateReefBiomass's
    # equivalent test -- would catch a wrong factor, wrong cutoff logic, or
    # missing dw weighting even if the factor-mirroring test above didn't.
    data(caribbean_3_model)
    params <- caribbean_3_model
    no_sp <- nrow(params@species_params)
    params@species_params$number_observed <- c(1000, NA, 500)[seq_len(no_sp)]
    params@species_params$number_cutoff <- rep(0, no_sp)
    observed_total <- sum(params@species_params$number_observed, na.rm = TRUE)

    result <- calibrateReefNumber(params)
    model_total <- 0
    for (sp_idx in which(!is.na(params@species_params$number_observed))) {
        model_total <- model_total + sum(result@initial_n[sp_idx, ] * result@dw)
    }
    expect_equal(model_total, observed_total)
})

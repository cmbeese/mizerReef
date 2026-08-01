test_that("calibrateReefBiomass rescales by exactly observed_total / model_total", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    cutoff <- params@species_params$biomass_cutoff
    observed <- params@species_params$biomass_observed
    observed_total <- sum(observed, na.rm = TRUE)

    model_total <- 0
    for (sp_idx in which(!is.na(observed))) {
        above_cutoff <- params@w >= cutoff[[sp_idx]]
        model_total <- model_total +
            sum((params@initial_n[sp_idx, ] * params@w * params@dw)[above_cutoff])
    }
    expected_factor <- observed_total / model_total

    expect_equal(calibrateReefBiomass(params), scaleReefModel(params, factor = expected_factor))
})

test_that("calibrateReefBiomass returns params unchanged when biomass_observed is absent", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@species_params$biomass_observed <- NULL
    expect_equal(calibrateReefBiomass(params), params)
})

test_that("calibrateReefBiomass returns params unchanged when biomass_observed is all NA", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@species_params$biomass_observed <- NA
    expect_equal(calibrateReefBiomass(params), params)
})

test_that("calibrateReefBiomass makes total model biomass (above cutoff) match total observed", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    observed_total <- sum(params@species_params$biomass_observed, na.rm = TRUE)
    cutoff <- params@species_params$biomass_cutoff

    result <- calibrateReefBiomass(params)
    model_total <- 0
    for (sp_idx in which(!is.na(params@species_params$biomass_observed))) {
        above_cutoff <- result@w >= cutoff[[sp_idx]]
        model_total <- model_total +
            sum((result@initial_n[sp_idx, ] * result@w * result@dw)[above_cutoff])
    }
    expect_equal(model_total, observed_total)
})

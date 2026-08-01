test_that("reducing grazer abundance increases tuned algae biomass, for fixed algae production", {
    # This is the core ecological behaviour the tuneUR()/tuneUR_cc() fix is
    # for: algal production on a reef is real primary production and is not
    # driven by grazer/herbivore demand, so cutting grazer abundance should
    # let algae accumulate (production unchanged, less gets eaten), not
    # silently reduce modelled algae production to match the new, lower
    # consumption (which is what the old tuneUR() did, by always resetting
    # algae_growth to the current total consumption).
    data(caribbean_3_model)
    params <- caribbean_3_model
    baseline <- suppressWarnings(tuneUR(params))
    baseline_growth <- baseline@other_params$algae$growth
    baseline_biomass <- algae_biomass(baseline)

    # "herbivores" is the only species with rho_algae > 0 in this model.
    fewer_grazers <- params
    fewer_grazers@initial_n["herbivores", ] <-
        fewer_grazers@initial_n["herbivores", ] * 0.5
    fewer_grazers <- suppressWarnings(tuneUR(fewer_grazers))

    # Algae production must be exactly the same fixed rate -- not retuned.
    expect_equal(
        fewer_grazers@other_params$algae$growth,
        baseline_growth
    )
    # But less consumption means more standing algae biomass.
    expect_gt(algae_biomass(fewer_grazers), baseline_biomass)
})

test_that("increasing grazer abundance decreases tuned algae biomass, for fixed algae production", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    baseline <- suppressWarnings(tuneUR(params))
    baseline_growth <- baseline@other_params$algae$growth
    baseline_biomass <- algae_biomass(baseline)

    more_grazers <- params
    more_grazers@initial_n["herbivores", ] <-
        more_grazers@initial_n["herbivores", ] * 2
    more_grazers <- suppressWarnings(tuneUR(more_grazers))

    expect_equal(
        more_grazers@other_params$algae$growth,
        baseline_growth
    )
    expect_lt(algae_biomass(more_grazers), baseline_biomass)
})

test_that("the same grazing-response pattern holds for tuneUR_cc (carrying-capacity variant)", {
    data(caribbean_3_model)
    cc_params <- setURcapacity(caribbean_3_model, cap = 1.5)
    baseline <- suppressWarnings(tuneUR_cc(cc_params))
    baseline_growth <- baseline@other_params$algae$growth
    baseline_biomass <- algae_biomass(baseline)

    fewer_grazers <- cc_params
    fewer_grazers@initial_n["herbivores", ] <-
        fewer_grazers@initial_n["herbivores", ] * 0.5
    fewer_grazers <- suppressWarnings(tuneUR_cc(fewer_grazers))

    expect_equal(
        fewer_grazers@other_params$algae$growth,
        baseline_growth
    )
    expect_gt(algae_biomass(fewer_grazers), baseline_biomass)
})

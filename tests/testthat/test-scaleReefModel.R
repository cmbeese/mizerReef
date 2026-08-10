test_that("scaleReefModel scales algae/detritus parameters as documented", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    old_algae_rho <- params@other_params$algae$rho
    old_algae_growth <- params@other_params$algae$growth
    old_detritus_rho <- params@other_params$detritus$rho
    old_detritus_external <- params@other_params$detritus$external

    result <- scaleReefModel(params, factor = 3)

    expect_equal(result@other_params$algae$rho, old_algae_rho / 3)
    expect_equal(result@other_params$detritus$rho, old_detritus_rho / 3)
    expect_equal(result@other_params$detritus$external, old_detritus_external * 3)
})

test_that("scaleReefModel leaves algae growth unscaled", {
    # algae_growth is a fixed, literature-informed production rate (see
    # getAlgaeProduction()) that the algae redesign deliberately holds
    # constant, independent of the model's abundance scale -- unlike rho,
    # which is an encounter-rate coefficient and must still be rescaled
    # to keep algae consumption pressure invariant under the rescale.
    data(caribbean_3_model)
    params <- caribbean_3_model
    old_algae_growth <- params@other_params$algae$growth

    result <- scaleReefModel(params, factor = 3)

    expect_equal(result@other_params$algae$growth, old_algae_growth)
})

test_that("scaleReefModel scales core mizer state (abundances, resource, search_vol) as mizer::scaleModel documents", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    old_n <- params@initial_n
    old_n_pp <- params@initial_n_pp
    old_search_vol <- params@search_vol
    old_cc_pp <- params@cc_pp

    result <- scaleReefModel(params, factor = 3)

    expect_equal(result@initial_n, old_n * 3)
    expect_equal(result@initial_n_pp, old_n_pp * 3)
    expect_equal(result@search_vol, old_search_vol / 3)
    expect_equal(result@cc_pp, old_cc_pp * 3)
})

test_that("scaleReefModel scales every unstructured resource's initial biomass, algae and detritus alike", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    old_algae_biomass <- algae_biomass(params)
    old_detritus_biomass <- detritus_biomass(params)

    result <- scaleReefModel(params, factor = 3)
    expect_equal(algae_biomass(result), old_algae_biomass * 3)
    expect_equal(detritus_biomass(result), old_detritus_biomass * 3)
})

test_that("scaleReefModel with factor = 1 leaves algae/detritus parameters unchanged", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- scaleReefModel(params, factor = 1)

    expect_equal(result@other_params$algae, params@other_params$algae)
    expect_equal(result@other_params$detritus, params@other_params$detritus)
})

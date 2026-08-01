test_that("rescale_algae multiplies biomass and divides rho/rho_algae by the same factor", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    old_biomass <- algae_biomass(params)
    old_rho <- params@other_params$algae$rho
    old_rho_algae <- params@species_params$rho_algae

    result <- rescale_algae(params, factor = 4)

    expect_equal(algae_biomass(result), old_biomass * 4)
    expect_equal(result@other_params$algae$rho, old_rho / 4)
    expect_equal(result@species_params$rho_algae, old_rho_algae / 4)
})

test_that("rescale_algae keeps total algae consumption (rate * biomass) unchanged", {
    # This is the documented purpose ("without changing anything else"):
    # biomass and mass-specific consumption rate move in opposite
    # directions so their product -- the total flux -- is invariant.
    data(caribbean_3_model)
    params <- caribbean_3_model
    old_total <- algae_consumption(params) * algae_biomass(params)

    result <- rescale_algae(params, factor = 3.7)
    new_total <- algae_consumption(result) * algae_biomass(result)

    expect_equal(new_total, old_total)
})

test_that("rescale_algae with factor = 1 is a no-op", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_equal(rescale_algae(params, 1), params)
})

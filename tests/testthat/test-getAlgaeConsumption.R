test_that("getAlgaeConsumption matches an independently computed per-species vector", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    feeding_level <- getFeedingLevel(params)
    rho <- params@other_params$algae$rho
    n <- params@initial_n
    B_A <- algae_biomass(params)

    expected <- vapply(seq_len(nrow(n)), function(sp) {
        sum(rho[sp, ] * n[sp, ] * (1 - feeding_level[sp, ]) * params@dw) * B_A
    }, numeric(1))
    names(expected) <- params@species_params$species

    expect_equal(as.numeric(getAlgaeConsumption(params)), as.numeric(expected))
})

test_that("getAlgaeConsumption agrees with the undiscounted rate for satiation = FALSE species", {
    # satiation = FALSE species have feeding_level forced to exactly 0 by
    # projectFeedingLevel.mizerReef() (unlimited intake), so (1 - f) = 1 and
    # getAlgaeConsumption()'s feeding-level-adjusted rate should collapse to
    # the plain rho*n*dw*B_A rate for exactly those species.
    data(caribbean_3_model)
    params <- caribbean_3_model
    unsatiated_idx <- which(!params@species_params$satiation)
    skip_if(length(unsatiated_idx) == 0, "no satiation = FALSE species in caribbean_3_model")

    rho <- params@other_params$algae$rho
    n <- params@initial_n
    B_A <- algae_biomass(params)
    raw <- vapply(unsatiated_idx, function(sp) {
        sum(rho[sp, ] * n[sp, ] * params@dw) * B_A
    }, numeric(1))

    expect_equal(as.numeric(getAlgaeConsumption(params))[unsatiated_idx], unname(raw))
})

test_that("getAlgaeConsumption returns one named entry per species", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_named(getAlgaeConsumption(params), params@species_params$species)
})

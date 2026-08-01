test_that("getDetritusConsumption matches an independently computed per-species vector", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    feeding_level <- getFeedingLevel(params)
    rho <- params@other_params$detritus$rho
    n <- params@initial_n
    B_D <- detritus_biomass(params)

    expected <- vapply(seq_len(nrow(n)), function(sp) {
        sum(rho[sp, ] * n[sp, ] * (1 - feeding_level[sp, ]) * params@dw) * B_D
    }, numeric(1))
    names(expected) <- params@species_params$species

    expect_equal(as.numeric(getDetritusConsumption(params)), as.numeric(expected))
})

test_that("getDetritusConsumption is zero for a species with no detritus interaction", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@other_params$detritus$rho[1, ] <- 0
    expect_equal(unname(getDetritusConsumption(params)[1]), 0)
})

test_that("getDetritusConsumption returns one named entry per species", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_named(getDetritusConsumption(params), params@species_params$species)
})

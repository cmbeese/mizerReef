test_that("algae_biomass returns exactly initial_n_other$algae", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_equal(algae_biomass(params), params@initial_n_other$algae)
})

test_that("algae_biomass reflects changes to initial_n_other$algae", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@initial_n_other$algae <- 12345
    expect_equal(algae_biomass(params), 12345)
})

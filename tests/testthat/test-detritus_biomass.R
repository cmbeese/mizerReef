test_that("detritus_biomass returns exactly initial_n_other$detritus", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_equal(detritus_biomass(params), params@initial_n_other$detritus)
})

test_that("detritus_biomass reflects changes to initial_n_other$detritus", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@initial_n_other$detritus <- 6789
    expect_equal(detritus_biomass(params), 6789)
})

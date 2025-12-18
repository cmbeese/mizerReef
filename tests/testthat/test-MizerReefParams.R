test_that("MizerReefParams returns correct class", {
    data(caribbean_3_model)
    sp <- caribbean_3_model@species_params
    params <- MizerReefParams(species_params = sp)
    expect_s4_class(params, "MizerReefParams")
})

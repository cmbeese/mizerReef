test_that("newReefParams returns correct class", {
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    data(tuning_profile)
    sp <- caribbean_3_species
    int <- caribbean_3_interaction
    result <- newReefParams(
        species_params = sp,
        interaction = int,
        method = "binned",
        method_params = tuning_profile
    )
    expect_s4_class(result, "mizerReef")
})

test_that("MizerReefParams issues deprecation warning", {
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    data(tuning_profile)
    expect_warning(
        params <- MizerReefParams(
            species_params = caribbean_3_species,
            interaction = caribbean_3_interaction,
            method = "binned",
            method_params = tuning_profile
        ),
        "deprecated"
    )
    expect_s4_class(params, "mizerReef")
})

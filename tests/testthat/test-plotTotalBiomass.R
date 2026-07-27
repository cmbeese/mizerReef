test_that("plotTotalBiomass(params) matches mizer::getBiomass() directly", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotTotalBiomass(params, return_data = TRUE)
    expected <- mizer::getBiomass(params)
    expect_equal(
        result$value[match(names(expected), result$Species)],
        unname(expected)
    )
})

test_that("plotTotalBiomass(sim) matches mizer::getBiomass() at the true last timestep", {
    # Regression test for the same off-by-one bug as plotTotalAbundance,
    # plus a crash: getBiomass.mizerReefSim() also returns algae/detritus
    # columns, which broke the (buggy) species-only logical subsetting.
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 5, progress_bar = FALSE)
    result <- plotTotalBiomass(sim, return_data = TRUE)

    expected <- mizer::getBiomass(sim)
    last_time <- as.character(max(as.numeric(dimnames(sim@n)$time)))
    fish_species <- caribbean_3_model@species_params$species
    expect_setequal(as.character(result$Species), fish_species)
    expect_equal(
        result$value[match(fish_species, result$Species)],
        unname(expected[last_time, fish_species])
    )
})

test_that("plotTotalBiomass respects an explicit non-contiguous species selection", {
    # Regression test for a fixed index-misalignment bug: selecting
    # species = c("predators", "inverts") (skipping "herbivores" in the
    # middle of the species order) used to attach herbivores' biomass to
    # a row labelled "inverts".
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotTotalBiomass(params, species = c("predators", "inverts"), return_data = TRUE)
    expected <- mizer::getBiomass(params)

    expect_setequal(as.character(result$Species), c("predators", "inverts"))
    expect_equal(result$value[result$Species == "predators"], unname(expected["predators"]))
    expect_equal(result$value[result$Species == "inverts"], unname(expected["inverts"]))
})

test_that("plotTotalBiomass returns a ggplot object", {
    data(caribbean_3_model)
    result <- plotTotalBiomass(caribbean_3_model)
    expect_s3_class(result, "ggplot")
})

test_that("plotlyTotalBiomass returns a plotly object", {
    data(caribbean_3_model)
    result <- plotlyTotalBiomass(caribbean_3_model)
    expect_s3_class(result, "plotly")
})

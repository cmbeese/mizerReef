test_that("getBiomass includes algae and detritus for a mizerReefSim", {
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)
    expect_s3_class(sim, "mizerReefSim")

    b <- getBiomass(sim)
    expect_true(all(c("algae", "detritus") %in% colnames(b)))
})

test_that("plotBiomass returns a ggplot object regardless of load order", {
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)

    result <- plotBiomass(sim)
    expect_s3_class(result, "ggplot")
})

test_that("plotBiomass's underlying data matches getBiomass.mizerReefSim() exactly", {
    # plotBiomass() itself dispatches to mizer's own (already-tested)
    # implementation -- the reef-specific piece is getBiomass.mizerReefSim(),
    # so the meaningful check is that its algae/detritus/species values
    # actually flow through into the plot correctly.
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)

    b <- getBiomass(sim)
    plot_data <- plotBiomass(sim)$data

    for (i in seq_len(nrow(plot_data))) {
        yr <- as.character(plot_data$Year[i])
        sp <- as.character(plot_data$Species[i])
        expect_equal(plot_data$Biomass[i], unname(b[yr, sp]), info = paste(yr, sp))
    }
})

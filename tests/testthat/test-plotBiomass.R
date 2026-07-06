test_that("getBiomass includes algae and detritus for a mizerReefSim", {
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)
    expect_s4_class(sim, "mizerReefSim")

    b <- getBiomass(sim)
    expect_true(all(c("algae", "detritus") %in% colnames(b)))
})

test_that("plotBiomass returns a ggplot object regardless of load order", {
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)

    result <- plotBiomass(sim)
    expect_s3_class(result, "ggplot")
})

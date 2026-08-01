test_that("plot2Productivity matches two separate plotProductivity() calls, rbind'd together", {
    data(caribbean_3_model)
    params1 <- caribbean_3_model
    sim2 <- project(caribbean_3_model, effort = 0.5, t_max = 2, progress_bar = FALSE)
    params2 <- sim2@params
    params2@initial_n <- sim2@n[dim(sim2@n)[1], , ]

    result <- plot2Productivity(params1, params2, name1 = "A", name2 = "B", return_data = TRUE)

    sf1 <- plotProductivity(params1, return_data = TRUE)
    sf2 <- plotProductivity(params2, return_data = TRUE)
    expect_equal(nrow(result), nrow(sf1) + nrow(sf2))
    expect_equal(result$value[result$Model == "A"], sf1$value)
    expect_equal(result$value[result$Model == "B"], sf2$value)
    expect_identical(levels(result$Model), c("A", "B"))
})

test_that("plot2Productivity works with MizerSim objects on both sides", {
    # Regression test for a fixed bug: plotProductivity()'s species-default
    # detection used missing(species), which was always FALSE once
    # plot2Productivity() forwarded species = species (even NULL),
    # breaking the MizerSim branch entirely ("object 'params' not found").
    data(caribbean_3_model)
    sim1 <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)
    sim2 <- project(caribbean_3_model, effort = 0.5, t_max = 2, progress_bar = FALSE)

    result <- plot2Productivity(sim1, sim2, name1 = "A", name2 = "B", return_data = TRUE)
    expect_setequal(as.character(result$Species), c("predators", "herbivores", "Total"))
})

test_that("plot2Productivity(include_inverts = TRUE) includes inverts", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plot2Productivity(params, params, include_inverts = TRUE, return_data = TRUE)
    expect_true("inverts" %in% as.character(result$Species))
})

test_that("plot2Productivity returns a ggplot object for both stack settings", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_s3_class(plot2Productivity(params, params, stack = FALSE), "ggplot")
    expect_s3_class(plot2Productivity(params, params, stack = TRUE), "ggplot")
})

test_that("plotly2Productivity returns a plotly object", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotly2Productivity(params, params)
    expect_s3_class(result, "plotly")
})

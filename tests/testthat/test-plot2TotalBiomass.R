test_that("plot2TotalBiomass matches two separate plotTotalBiomass() calls, rbind'd together", {
    data(caribbean_3_model)
    params1 <- caribbean_3_model
    sim2 <- project(caribbean_3_model, effort = 0.5, t_max = 2, progress_bar = FALSE)
    params2 <- sim2@params
    params2@initial_n <- sim2@n[dim(sim2@n)[1], , ]

    result <- plot2TotalBiomass(params1, params2, name1 = "A", name2 = "B", return_data = TRUE)

    sf1 <- plotTotalBiomass(params1, return_data = TRUE)
    sf2 <- plotTotalBiomass(params2, return_data = TRUE)
    expect_equal(result$value[result$Model == "A"], sf1$value)
    expect_equal(result$value[result$Model == "B"], sf2$value)
})

test_that("plot2TotalBiomass matches two direct plotTotalBiomass(sim) calls", {
    data(caribbean_3_model)
    sim1 <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)
    sim2 <- project(caribbean_3_model, effort = 0.5, t_max = 2, progress_bar = FALSE)

    result <- plot2TotalBiomass(sim1, sim2, name1 = "A", name2 = "B", return_data = TRUE)
    sf1 <- plotTotalBiomass(sim1, return_data = TRUE)
    sf2 <- plotTotalBiomass(sim2, return_data = TRUE)
    expect_equal(result$value[result$Model == "A"], sf1$value)
    expect_equal(result$value[result$Model == "B"], sf2$value)
})

test_that("plot2TotalBiomass returns a ggplot object for both stack settings", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_s3_class(plot2TotalBiomass(params, params, stack = FALSE), "ggplot")
    expect_s3_class(plot2TotalBiomass(params, params, stack = TRUE), "ggplot")
})

test_that("plotly2TotalBiomass returns a plotly object", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotly2TotalBiomass(params, params)
    expect_s3_class(result, "plotly")
})

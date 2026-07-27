test_that("plotSpectraPercentChange returns a ggplot object", {
    data(caribbean_3_model)
    sim1 <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)
    sim2 <- project(caribbean_3_model, effort = 0.5, t_max = 2,
                    progress_bar = FALSE)

    result <- plotSpectraPercentChange(sim1, sim2, power = 1)
    expect_s3_class(result, "ggplot")
})

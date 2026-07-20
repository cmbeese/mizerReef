test_that("plotRefugeDensity returns a ggplot object", {
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)
    result <- plotRefugeDensity(sim, return_data = FALSE)
    expect_s3_class(result, "ggplot")
})

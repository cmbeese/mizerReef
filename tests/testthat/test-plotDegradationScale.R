test_that("plotDegradationScale returns a ggplot object", {
    data(caribbean_3_model)
    # Set up degradation parameters (deg_scale is a required matrix)
    deg_scale <- matrix(c(0.8, 0.85, 0.9), nrow = 1, ncol = 3)
    params <- setDegradation(
        caribbean_3_model,
        deg_scale = deg_scale,
        degrade = TRUE
    )
    result <- plotDegradationScale(object = params)
    expect_s3_class(result, "ggplot")
})

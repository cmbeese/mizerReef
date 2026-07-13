test_that("plotDegradationScale returns a ggplot object", {
    data(caribbean_3_model)
    data(rubble_scale)
    params <- setDegradation(
        caribbean_3_model,
        deg_scale = rubble_scale,
        degrade = TRUE
    )
    result <- plotDegradationScale(object = params)
    expect_s3_class(result, "ggplot")
})

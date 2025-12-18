test_that("plotDegradationScale returns a ggplot object", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotDegradationScale(object = params)
    expect_s3_class(result, "ggplot")
})

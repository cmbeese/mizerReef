test_that("plotDegScale returns a ggplot object", {
    result <- plotDegScale()
    expect_s3_class(result, "ggplot")
})

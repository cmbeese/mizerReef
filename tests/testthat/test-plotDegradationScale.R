test_that("plotDegradationScale with a named built-in trajectory matches independently melting that data object", {
    data(rubble_scale)
    result <- plotDegradationScale(trajectory = "rubble", return_data = TRUE)

    expected <- reshape2::melt(rubble_scale)
    colnames(expected) <- c("SizeBin", "Year", "Scaling")
    expect_equal(as.character(result$SizeBin), as.character(expected$SizeBin))
    expect_equal(result$Year, expected$Year)
    expect_equal(result$Scaling, expected$Scaling)
})

test_that("plotDegradationScale with a custom matrix trajectory matches independently melting it", {
    custom_traj <- matrix(c(0.5, 0.6, 0.7, 0.8),
        nrow = 2, dimnames = list(c("bin1", "bin2"), c("1", "2"))
    )
    result <- plotDegradationScale(trajectory = custom_traj, return_data = TRUE)

    expected <- reshape2::melt(custom_traj)
    expect_equal(as.character(result$SizeBin), as.character(expected$Var1))
    expect_equal(result$Scaling, expected$value)
})

test_that("plotDegradationScale with a MizerParams object pulls its stored deg_scale", {
    data(caribbean_3_model)
    data(rubble_scale)
    params <- setDegradation(caribbean_3_model, deg_scale = rubble_scale, degrade = TRUE)

    result <- plotDegradationScale(object = params, return_data = TRUE)
    expected <- reshape2::melt(rubble_scale)
    expect_equal(result$Scaling, expected[, 3])
})

test_that("plotDegradationScale errors when neither object nor trajectory is provided", {
    expect_error(plotDegradationScale(), "Either 'object' or 'trajectory' must be provided")
})

test_that("plotDegradationScale errors on an invalid trajectory name", {
    expect_error(plotDegradationScale(trajectory = "bogus"), "Invalid trajectory")
})

test_that("plotDegradationScale errors when the object has no degradation scaling configured", {
    data(caribbean_3_model)
    expect_error(
        plotDegradationScale(object = caribbean_3_model),
        "No degradation scaling or trajectory found"
    )
})

test_that("plotDegradationScale returns a ggplot object when return_data = FALSE", {
    data(caribbean_3_model)
    data(rubble_scale)
    params <- setDegradation(caribbean_3_model, deg_scale = rubble_scale, degrade = TRUE)
    result <- plotDegradationScale(object = params)
    expect_s3_class(result, "ggplot")
})

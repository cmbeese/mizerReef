test_that("plotDegScale's data matches independently melting each built-in trajectory", {
    data(rubble_scale)
    data(algae_scale)
    data(recovery_scale)
    trajectories <- list(rubble = rubble_scale, algae = algae_scale, recovery = recovery_scale)

    result <- plotDegScale(return_data = TRUE)
    expect_setequal(as.character(unique(result$Trajectory)), names(trajectories))

    for (traj in names(trajectories)) {
        expected <- reshape2::melt(trajectories[[traj]])
        colnames(expected) <- c("SizeBin", "Year", "Scaling")
        expected$SizeBin <- as.character(expected$SizeBin)
        expected <- expected[order(expected$SizeBin, expected$Year), ]
        rownames(expected) <- NULL

        actual <- result[result$Trajectory == traj, c("SizeBin", "Year", "Scaling")]
        actual$SizeBin <- as.character(actual$SizeBin)
        actual <- actual[order(actual$SizeBin, actual$Year), ]
        rownames(actual) <- NULL

        expect_equal(actual, expected, info = traj)
    }
})

test_that("plotDegScale returns a ggplot object when return_data = FALSE", {
    result <- plotDegScale()
    expect_s3_class(result, "ggplot")
})

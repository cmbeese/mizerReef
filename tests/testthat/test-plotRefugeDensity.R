test_that("plotRefugeDensity's data matches getDegrade() reshaped to long format", {
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 3, progress_bar = FALSE)

    result <- plotRefugeDensity(sim, return_data = TRUE)

    degrade <- getDegrade(sim)
    expected <- as.data.frame(as.table(degrade))
    colnames(expected) <- c("time", "size_bin", "refuge_density")
    expected$time <- as.numeric(as.character(expected$time))
    expected$size_bin <- as.numeric(as.character(expected$size_bin))

    expect_equal(result, expected)
})

test_that("plotRefugeDensity returns a ggplot object when return_data = FALSE", {
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)
    result <- plotRefugeDensity(sim, return_data = FALSE)
    expect_s3_class(result, "ggplot")
})

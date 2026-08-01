test_that("plotAlgaeConsumption's underlying data matches getAlgaeConsumption() filtered to >1% of total", {
    # plotAlgaeConsumption() has no return_data argument, but the ggplot
    # object's underlying data frame is accessible via result$data.
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotAlgaeConsumption(params)

    consumption <- getAlgaeConsumption(params)
    total <- sum(consumption)
    expected <- consumption[consumption > total / 100]

    plot_data <- result$data
    expect_setequal(as.character(plot_data$Consumer), names(expected))
    expect_equal(
        plot_data$Rate[match(names(expected), plot_data$Consumer)],
        unname(expected)
    )
})

test_that("plotAlgaeConsumption returns a valid ggplot when no species consumes algae", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@other_params$algae$rho[] <- 0

    result <- plotAlgaeConsumption(params)
    expect_s3_class(result, "ggplot")
    expect_equal(nrow(result$data), 0)
})

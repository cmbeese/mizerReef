test_that("plotDetritusConsumption's underlying data matches getDetritusConsumption() filtered to >1% of total", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotDetritusConsumption(params)

    consumption <- getDetritusConsumption(params)
    total <- sum(consumption)
    expected <- consumption[consumption > total / 100]

    plot_data <- result$data
    expect_setequal(as.character(plot_data$Consumer), names(expected))
    expect_equal(
        plot_data$Rate[match(names(expected), plot_data$Consumer)],
        unname(expected)
    )
})

test_that("plotDetritusConsumption shows the 'no consumers' message when no species consumes detritus", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@other_params$detritus$rho[] <- 0

    result <- plotDetritusConsumption(params)
    expect_s3_class(result, "ggplot")
    expect_equal(result$labels$title, "No detritus consumers with >1% of total consumption")
})

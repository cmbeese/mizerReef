test_that("plotDetritusProduction's underlying data matches getDetritusProduction()", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotDetritusProduction(params)

    expected <- getDetritusProduction(params)
    plot_data <- result$data

    expect_setequal(as.character(plot_data$Source), names(expected))
    expect_equal(
        plot_data$Rate[match(names(expected), plot_data$Source)],
        unname(expected)
    )
})

test_that("plotTotalAbundance(params) matches mizer::getN() directly", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotTotalAbundance(params, return_data = TRUE)
    expected <- mizer::getN(params)
    expect_equal(
        result$value[match(names(expected), result$Species)],
        unname(expected)
    )
})

test_that("plotTotalAbundance(sim) matches mizer::getN() at the true last timestep", {
    # Regression test for a fixed off-by-one bug: end_time (a time value)
    # was used directly as a row position, silently returning the
    # second-to-last saved timestep instead of the true last one.
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 5, progress_bar = FALSE)
    result <- plotTotalAbundance(sim, return_data = TRUE)

    expected <- mizer::getN(sim)
    last_time <- as.character(max(as.numeric(dimnames(sim@n)$time)))
    expect_equal(
        result$value[match(names(expected[last_time, ]), result$Species)],
        unname(expected[last_time, ])
    )
})

test_that("plotTotalAbundance respects an explicit non-contiguous species selection", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotTotalAbundance(params, species = c("predators", "inverts"), return_data = TRUE)
    expected <- mizer::getN(params)

    expect_setequal(as.character(result$Species), c("predators", "inverts"))
    expect_equal(
        result$value[match(names(expected)[c(1, 3)], result$Species)],
        unname(expected[c(1, 3)])
    )
})

test_that("plotTotalAbundance returns a ggplot object for both stack settings", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_s3_class(plotTotalAbundance(params, stack = FALSE), "ggplot")
    expect_s3_class(plotTotalAbundance(params, stack = TRUE), "ggplot")
})

test_that("plotlyTotalAbundance returns a plotly object", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotlyTotalAbundance(params)
    expect_s3_class(result, "plotly")
})

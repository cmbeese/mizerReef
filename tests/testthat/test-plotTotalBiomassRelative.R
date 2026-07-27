test_that("plotTotalBiomassRelative percent_change matches 100*(new-old)/old computed from getBiomass() directly", {
    # Regression test for a fixed bug: percent_change used to omit the *100,
    # plotting a raw fraction on an axis labelled "%".
    data(caribbean_3_model)
    params1 <- caribbean_3_model
    sim2 <- project(caribbean_3_model, effort = 0.5, t_max = 3, progress_bar = FALSE)
    params2 <- sim2@params
    params2@initial_n <- sim2@n[dim(sim2@n)[1], , ]

    result <- plotTotalBiomassRelative(params1, params2, diff_method = "percent_change", return_data = TRUE)
    b_a <- mizer::getBiomass(params1)
    b_b <- mizer::getBiomass(params2)
    expected <- 100 * (b_b - b_a) / b_a

    expect_equal(
        result$rel_diff[match(as.character(result$Species), names(expected))],
        unname(expected[as.character(result$Species)])
    )
})

test_that("plotTotalBiomassRelative rel_diff matches (new-old)/(old+new) computed from getBiomass() directly", {
    data(caribbean_3_model)
    params1 <- caribbean_3_model
    sim2 <- project(caribbean_3_model, effort = 0.5, t_max = 3, progress_bar = FALSE)
    params2 <- sim2@params
    params2@initial_n <- sim2@n[dim(sim2@n)[1], , ]

    result <- plotTotalBiomassRelative(params1, params2, diff_method = "rel_diff", return_data = TRUE)
    b_a <- mizer::getBiomass(params1)
    b_b <- mizer::getBiomass(params2)
    expected <- (b_b - b_a) / (b_a + b_b)

    expect_equal(
        result$rel_diff[match(as.character(result$Species), names(expected))],
        unname(expected[as.character(result$Species)])
    )
})

test_that("plotTotalBiomassRelative errors on an invalid diff_method", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(
        plotTotalBiomassRelative(params, params, diff_method = "bogus"),
        "diff_method should be either"
    )
})

test_that("plotTotalBiomassRelative returns a ggplot object", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotTotalBiomassRelative(params, params, diff_method = "percent_change")
    expect_s3_class(result, "ggplot")
})

test_that("plotlyTotalBiomassRelative returns a plotly object", {
    # Regression test for a fixed bug: this wrapper called the non-existent
    # function "TotalBiomassRelative" instead of "plotTotalBiomassRelative".
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotlyTotalBiomassRelative(params, params, diff_method = "percent_change")
    expect_s3_class(result, "plotly")
})

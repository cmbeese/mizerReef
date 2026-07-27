test_that("plotProductivityRelative percent_change matches 100*(new-old)/old computed from getProductivity() directly", {
    # Regression test for a fixed bug: percent_change used to omit the *100,
    # plotting a raw fraction on an axis labelled "%".
    data(caribbean_3_model)
    params1 <- caribbean_3_model
    sim2 <- project(caribbean_3_model, effort = 0.5, t_max = 3, progress_bar = FALSE)
    params2 <- sim2@params
    params2@initial_n <- sim2@n[dim(sim2@n)[1], , ]

    result <- plotProductivityRelative(params1, params2, diff_method = "percent_change", return_data = TRUE)
    p_a <- getProductivity(params1)
    p_b <- getProductivity(params2)
    expected <- 100 * (p_b - p_a) / p_a

    expect_equal(
        result$rel_diff[match(as.character(result$Species), names(expected))],
        unname(expected[as.character(result$Species)])
    )
})

test_that("plotProductivityRelative rel_diff matches (new-old)/(old+new) computed from getProductivity() directly", {
    data(caribbean_3_model)
    params1 <- caribbean_3_model
    sim2 <- project(caribbean_3_model, effort = 0.5, t_max = 3, progress_bar = FALSE)
    params2 <- sim2@params
    params2@initial_n <- sim2@n[dim(sim2@n)[1], , ]

    result <- plotProductivityRelative(params1, params2, diff_method = "rel_diff", return_data = TRUE)
    p_a <- getProductivity(params1)
    p_b <- getProductivity(params2)
    expected <- (p_b - p_a) / (p_a + p_b)

    expect_equal(
        result$rel_diff[match(as.character(result$Species), names(expected))],
        unname(expected[as.character(result$Species)])
    )
})

test_that("plotProductivityRelative errors on an invalid diff_method", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(
        plotProductivityRelative(params, params, diff_method = "bogus"),
        "diff_method should be either"
    )
})

test_that("plotProductivityRelative returns a ggplot object", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotProductivityRelative(params, params, diff_method = "percent_change")
    expect_s3_class(result, "ggplot")
})

test_that("plotlyProductivityRelative returns a plotly object", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotlyProductivityRelative(params, params, diff_method = "percent_change")
    expect_s3_class(result, "plotly")
})

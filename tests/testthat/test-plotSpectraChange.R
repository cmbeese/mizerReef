test_that("plotSpectraChange computes percent change matching an independent calculation from mizer::plotSpectra()", {
    # rel_diff = 100 * (N2(w) - N1(w)) / N1(w), hand-computed here from
    # mizer::plotSpectra(return_data = TRUE) called directly on each
    # simulation, not by re-deriving plotSpectraChange()'s own code.
    data(caribbean_3_model)
    sim1 <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)
    sim2 <- project(caribbean_3_model, effort = 0.5, t_max = 2, progress_bar = FALSE)

    sf1 <- mizer::plotSpectra(sim1, power = 1, return_data = TRUE)
    sf2 <- mizer::plotSpectra(sim2, power = 1, return_data = TRUE)
    names(sf1)[2] <- "value"
    names(sf2)[2] <- "value"
    sf <- dplyr::left_join(sf1, sf2, by = c("w", "Legend"))
    expected <- 100 * (sf$value.y - sf$value.x) / sf$value.x

    result <- plotSpectraChange(sim1, sim2, power = 1, return_data = TRUE)
    expect_equal(result$rel_diff, expected)
})

test_that("plotSpectraChange with use_percent = FALSE gives the raw relative proportion, exactly 1/100 of the percent value", {
    data(caribbean_3_model)
    sim1 <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)
    sim2 <- project(caribbean_3_model, effort = 0.5, t_max = 2, progress_bar = FALSE)

    pct <- plotSpectraChange(sim1, sim2, power = 1, use_percent = TRUE, return_data = TRUE)
    prop <- plotSpectraChange(sim1, sim2, power = 1, use_percent = FALSE, return_data = TRUE)

    expect_equal(pct$rel_diff, 100 * prop$rel_diff)
})

test_that("plotSpectraChange uses use_percent = TRUE by default", {
    data(caribbean_3_model)
    sim1 <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)
    sim2 <- project(caribbean_3_model, effort = 0.5, t_max = 2, progress_bar = FALSE)

    default_data <- plotSpectraChange(sim1, sim2, power = 1, return_data = TRUE)
    explicit_pct_data <- plotSpectraChange(sim1, sim2, power = 1, use_percent = TRUE, return_data = TRUE)
    expect_equal(default_data, explicit_pct_data)
})

test_that("plotSpectraChange labels the y-axis according to use_percent", {
    data(caribbean_3_model)
    sim1 <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)
    sim2 <- project(caribbean_3_model, effort = 0.5, t_max = 2, progress_bar = FALSE)

    pct_plot <- plotSpectraChange(sim1, sim2, power = 1, use_percent = TRUE)
    prop_plot <- plotSpectraChange(sim1, sim2, power = 1, use_percent = FALSE)

    expect_s3_class(pct_plot, "ggplot")
    expect_equal(pct_plot$labels$y, "% Change in Biomass")
    expect_equal(prop_plot$labels$y, "Relative Change in Biomass")
})

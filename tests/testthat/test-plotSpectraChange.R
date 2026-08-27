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
    expect_equal(pct_plot$labels$y, "% Change in Biomass density")
    expect_equal(prop_plot$labels$y, "Relative Change in Biomass density")
})

test_that("plotSpectraChange accepts mizer 3.3's biomass/per_log_size in place of power", {
    # mizer 3.3 replaced `power` with the two independent flags `biomass` and
    # `per_log_size`, with power = biomass + per_log_size. Passing the flags
    # must give the same numbers as the equivalent power, and the y-axis label
    # must name the quantity actually plotted rather than assuming biomass.
    data(caribbean_3_model)
    sim1 <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)
    sim2 <- project(caribbean_3_model, effort = 0.5, t_max = 2, progress_bar = FALSE)

    expect_equal(
        plotSpectraChange(sim1, sim2, biomass = TRUE, return_data = TRUE),
        plotSpectraChange(sim1, sim2, power = 1, return_data = TRUE)
    )
    expect_equal(
        plotSpectraChange(sim1, sim2, biomass = TRUE, per_log_size = TRUE,
                          return_data = TRUE),
        plotSpectraChange(sim1, sim2, power = 2, return_data = TRUE)
    )

    number_plot <- plotSpectraChange(sim1, sim2, biomass = FALSE)
    expect_equal(number_plot$labels$y, "% Change in Number density")
})

test_that("plotSpectraChange defaults to the biomass density, as plotSpectra() does", {
    data(caribbean_3_model)
    sim1 <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)
    sim2 <- project(caribbean_3_model, effort = 0.5, t_max = 2, progress_bar = FALSE)

    expect_equal(
        plotSpectraChange(sim1, sim2, return_data = TRUE),
        plotSpectraChange(sim1, sim2, power = 1, return_data = TRUE)
    )
})

test_that("plotSpectraChange follows plotSpectra() onto a length axis", {
    # `...` goes to plotSpectra(), so size_axis = "l" reaches it and the
    # returned data frame has an `l` column instead of `w`. The join, the
    # x aesthetic and the axis label all have to follow.
    data(caribbean_3_model)
    d <- plotSpectraChange(caribbean_3_model, caribbean_3_model,
                           size_axis = "l", return_data = TRUE)
    expect_true("l" %in% names(d))
    expect_false("w" %in% names(d))

    p <- plotSpectraChange(caribbean_3_model, caribbean_3_model,
                           size_axis = "l")
    expect_s3_class(p, "ggplot")
    expect_equal(p$labels$x, "Length [cm]")

    expect_equal(
        plotSpectraChange(caribbean_3_model, caribbean_3_model)$labels$x,
        "Weight [g]"
    )
})

test_that("spectra_quantity() names the quantity exactly as mizer's own column does", {
    # plotSpectraChange() reads the plotted quantity off plotSpectra()'s value
    # column, and falls back to spectra_quantity() when an extension package's
    # method renamed that column (mizerMR calls it "value"). The fallback has
    # to agree with mizer, so check it against mizer over every combination of
    # `power`, `biomass` and `per_log_size`, units stripped.
    combos <- list(
        list(), list(power = 0), list(power = 1), list(power = 2),
        list(power = 1.5), list(biomass = TRUE), list(biomass = FALSE),
        list(biomass = FALSE, per_log_size = TRUE),
        list(biomass = TRUE, per_log_size = TRUE),
        list(per_log_size = TRUE), list(power = 1, biomass = FALSE)
    )
    for (args in combos) {
        from_mizer <- sub(
            " \\[.*\\]$", "",
            names(do.call(mizer::plotSpectra,
                          c(list(mizer::NS_params, return_data = TRUE),
                            args)))[2]
        )
        expect_equal(do.call(spectra_quantity, args), from_mizer,
                     info = paste(names(args), unlist(args),
                                  sep = "=", collapse = ", "))
    }
})

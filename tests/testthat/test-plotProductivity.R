test_that("plotProductivity(params) matches getProductivity() and excludes inverts by default", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotProductivity(params, return_data = TRUE)
    expected <- getProductivity(params)

    expect_setequal(as.character(result$Species), c("predators", "herbivores"))
    expect_equal(
        result$value[match(c("predators", "herbivores"), result$Species)],
        unname(expected[c("predators", "herbivores")])
    )
})

test_that("plotProductivity(params, include_inverts = TRUE) includes all species", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotProductivity(params, include_inverts = TRUE, return_data = TRUE)
    expected <- getProductivity(params)

    expect_setequal(as.character(result$Species), names(expected))
    expect_equal(
        result$value[match(names(expected), result$Species)],
        unname(expected)
    )
})

test_that("plotProductivity(params) with an explicit species argument overrides include_inverts", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotProductivity(params, species = c("predators", "inverts"), return_data = TRUE)
    expect_setequal(as.character(result$Species), c("predators", "inverts"))
})

test_that("plotProductivity(sim) matches getProductivity() over the full time range, excluding inverts by default", {
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 3, progress_bar = FALSE)
    result <- plotProductivity(sim, return_data = TRUE)

    p <- getProductivity(sim, time_range = 0:3)
    p <- cbind(p, Total = rowSums(p))
    expected <- reshape2::melt(p)
    names(expected) <- c("Year", "Species", "Productivity")
    expected <- expected[expected$Species %in% c("Total", "predators", "herbivores"), ]

    expect_equal(nrow(result), nrow(expected))
    for (yr in unique(expected$Year)) {
        for (sp in unique(as.character(expected$Species))) {
            actual_val <- result$Productivity[result$Year == yr & result$Species == sp]
            expected_val <- expected$Productivity[expected$Year == yr & expected$Species == sp]
            expect_equal(actual_val, expected_val, info = paste(yr, sp))
        }
    }
})

test_that("plotProductivity(sim, include_inverts = TRUE) includes inverts", {
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)
    result <- plotProductivity(sim, include_inverts = TRUE, return_data = TRUE)
    expect_true("inverts" %in% as.character(result$Species))
})

test_that("plotProductivity(sim) errors when start_time >= end_time", {
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 2, progress_bar = FALSE)
    expect_error(
        plotProductivity(sim, start_time = 2, end_time = 0, return_data = TRUE),
        "start_time must be less than end_time"
    )
})

test_that("plotProductivity returns a ggplot object", {
    data(caribbean_3_model)
    result <- plotProductivity(caribbean_3_model)
    expect_s3_class(result, "ggplot")
})

test_that("plotlyProductivity returns a plotly object", {
    data(caribbean_3_model)
    result <- plotlyProductivity(caribbean_3_model)
    expect_s3_class(result, "plotly")
})

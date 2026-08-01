test_that("plotRelativeContribution's rel column is each species' share of the metric total, and excludes inverts by default", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotRelativeContribution(params, return_data = TRUE)

    expect_setequal(as.character(result$Species), c("predators", "herbivores"))
    expect_setequal(as.character(result$Metric), c("Abundance", "Biomass", "Productivity"))

    abd <- mizer::getN(params)[c("predators", "herbivores")]
    biom <- mizer::getBiomass(params)[c("predators", "herbivores")]
    prod <- getProductivity(params)[c("predators", "herbivores")]

    for (m in c("Abundance", "Biomass", "Productivity")) {
        sub <- result[result$Metric == m, ]
        expected_value <- switch(m, Abundance = abd, Biomass = biom, Productivity = prod)
        expect_equal(
            sub$value[match(names(expected_value), sub$Species)],
            unname(expected_value)
        )
        expect_equal(sum(sub$rel), 1)
        expect_equal(
            sub$rel[match(names(expected_value), sub$Species)],
            unname(expected_value) / sum(expected_value)
        )
    }
})

test_that("plotRelativeContribution(include_inverts = TRUE) includes inverts in all three metrics", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotRelativeContribution(params, include_inverts = TRUE, return_data = TRUE)

    expect_setequal(as.character(result$Species), c("predators", "herbivores", "inverts"))
    for (m in unique(as.character(result$Metric))) {
        expect_equal(sum(result$rel[result$Metric == m]), 1)
    }
})

test_that("plotRelativeContribution returns a ggplot object", {
    data(caribbean_3_model)
    result <- plotRelativeContribution(caribbean_3_model)
    expect_s3_class(result, "ggplot")
})

test_that("plotlyRelativeContribution returns a plotly object", {
    data(caribbean_3_model)
    result <- plotlyRelativeContribution(caribbean_3_model)
    expect_s3_class(result, "plotly")
})

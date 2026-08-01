test_that("plotRefugeProfile equals 1 - getVulnerable() for species with refuge_user = TRUE", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    refuge <- 1 - getVulnerable(params)
    refuge_users <- params@species_params$species[params@species_params$refuge_user]

    result <- plotRefugeProfile(params, return_data = TRUE, all.sizes = TRUE)
    expect_setequal(as.character(result$Species), refuge_users)
    for (sp in refuge_users) {
        expect_equal(
            result$value[result$Species == sp],
            unname(refuge[sp, ]),
            info = sp
        )
    }
})

test_that("plotRefugeProfile excludes species with refuge_user = FALSE", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    non_users <- params@species_params$species[!params@species_params$refuge_user]
    expect_true(length(non_users) > 0)

    result <- plotRefugeProfile(params, return_data = TRUE)
    expect_true(!any(non_users %in% as.character(result$Species)))
})

test_that("plotRefugeProfile's length bins use each species' own a/b conversion", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotRefugeProfile(params, return_data = TRUE, all.sizes = TRUE)

    for (sp in unique(as.character(result$Species))) {
        a <- params@species_params[sp, "a"]
        b <- params@species_params[sp, "b"]
        w <- result$w[result$Species == sp]
        expected_l <- (w / a)^(1 / b)
        expect_equal(result$l[result$Species == sp], expected_l, info = sp)
    }
})

test_that("plotRefugeProfile works on a MizerSim object (using its steady state) and warns", {
    # Regression test for a fixed bug: an unconditional params <- object
    # right after the if/else block undid the MizerSim branch's correct
    # params <- object@params extraction, crashing with
    # "no slot of name 'species_params'".
    data(caribbean_3_model)
    params <- caribbean_3_model
    sim <- project(params, t_max = 2, progress_bar = FALSE)

    expect_warning(
        result <- plotRefugeProfile(sim, return_data = TRUE, all.sizes = TRUE),
        "steady state"
    )
    expected <- plotRefugeProfile(params, return_data = TRUE, all.sizes = TRUE)
    expect_equal(result, expected)
})

test_that("plotRefugeProfile returns a ggplot object", {
    data(caribbean_3_model)
    result <- plotRefugeProfile(caribbean_3_model)
    expect_s3_class(result, "ggplot")
})

test_that("plotlyRefugeProfile returns a plotly object", {
    data(caribbean_3_model)
    result <- plotlyRefugeProfile(caribbean_3_model)
    expect_s3_class(result, "plotly")
})

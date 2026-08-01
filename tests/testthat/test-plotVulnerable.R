test_that("plotVulnerable(params) matches getVulnerable() for species with refuge_user = TRUE", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    vul <- getVulnerable(params)
    refuge_users <- params@species_params$species[params@species_params$refuge_user]

    result <- plotVulnerable(params, return_data = TRUE, all.sizes = TRUE)
    expect_setequal(as.character(result$Species), refuge_users)
    for (sp in refuge_users) {
        expect_equal(
            result$value[result$Species == sp],
            unname(vul[sp, ]),
            info = sp
        )
    }
})

test_that("plotVulnerable excludes species with refuge_user = FALSE", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    non_users <- params@species_params$species[!params@species_params$refuge_user]
    expect_true(length(non_users) > 0)

    result <- plotVulnerable(params, return_data = TRUE)
    expect_true(!any(non_users %in% as.character(result$Species)))
})

test_that("plotVulnerable(sim) matches getVulnerable(sim) at the true last timestep by default", {
    # Regression test for a fixed bug: t <- time_step (NULL by default)
    # unconditionally overwrote the computed last-timestep default,
    # crashing plotVulnerable(sim) whenever time_step wasn't given.
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 3, progress_bar = FALSE)
    last_t <- max(as.numeric(dimnames(sim@n)$time))
    vul <- getVulnerable(sim, time_range = last_t)

    result <- plotVulnerable(sim, return_data = TRUE, all.sizes = TRUE)
    expect_equal(
        result$value[result$Species == "herbivores"],
        unname(vul["herbivores", ])
    )
})

test_that("plotVulnerable(sim) with an explicit time_step matches getVulnerable() at that time", {
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 3, progress_bar = FALSE)
    vul_t1 <- getVulnerable(sim, time_range = 1)

    result <- plotVulnerable(sim, time_step = 1, return_data = TRUE, all.sizes = TRUE)
    expect_equal(
        result$value[result$Species == "herbivores"],
        unname(vul_t1["herbivores", ])
    )
})

test_that("plotVulnerable with all.sizes = FALSE restricts to each species' own size range", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- plotVulnerable(params, return_data = TRUE, all.sizes = FALSE)

    for (sp in unique(as.character(result$Species))) {
        w_range <- range(result$w[result$Species == sp])
        expect_gte(w_range[1], params@species_params[sp, "w_min"])
        expect_lte(w_range[2], params@species_params[sp, "w_max"])
    }
})

test_that("plotVulnerable returns a ggplot object", {
    data(caribbean_3_model)
    result <- plotVulnerable(caribbean_3_model)
    expect_s3_class(result, "ggplot")
})

test_that("plotlyVulnerable returns a plotly object", {
    data(caribbean_3_model)
    result <- plotlyVulnerable(caribbean_3_model)
    expect_s3_class(result, "plotly")
})

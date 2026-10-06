# mizer's own getDiet() method, which ignores predation refuge
mizer_getDiet <- utils::getS3method("getDiet", "MizerParams")

# For each fish prey, the biomass all predators consume and the biomass that
# prey loses to predation. These must agree.
prey_balance <- function(params, diet) {
    n <- params@initial_n
    pred_mort <- getPredMort(params)
    dw <- rep(params@dw, each = nrow(n))
    sapply(params@species_params$species, function(prey) c(
        consumed = sum(diet[, , prey] * n * dw),
        lost = sum(pred_mort[prey, ] * n[prey, ] * params@w * params@dw)
    ))
}

# Total consumption equals encounter x (1 - feeding level), wherever fish exist
expected_total <- function(params) {
    total <- getEncounter(params) * (1 - getFeedingLevel(params))
    total[params@initial_n <= 0] <- 0
    total
}

test_that("getDiet() consumption equals predation losses in both bundled models", {
    data(caribbean_3_model, caribbean_10_model)
    for (params in list(caribbean_3_model, caribbean_10_model)) {
        bal <- prey_balance(params, getDiet(params, proportion = FALSE))
        expect_equal(bal["consumed", ], bal["lost", ], tolerance = 1e-8)
    }
})

test_that("mizer's own getDiet() method overstates prey that use refuge", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    bal <- prey_balance(params, mizer_getDiet(params, proportion = FALSE))
    # predators use refuge, inverts don't
    expect_gt(bal["consumed", "predators"], 1.1 * bal["lost", "predators"])
    expect_equal(bal["consumed", "inverts"], bal["lost", "inverts"], tolerance = 1e-8)
})

test_that("getDiet() total equals encounter x (1 - feeding level)", {
    data(caribbean_3_model, caribbean_10_model)
    for (params in list(caribbean_3_model, caribbean_10_model)) {
        total <- rowSums(getDiet(params, proportion = FALSE), dims = 2)
        expect_equal(total, expected_total(params), tolerance = 1e-8,
                     ignore_attr = TRUE)
    }
})

test_that("getDiet() leaves the diet of predators not blocked by refuge unchanged", {
    data(caribbean_10_model)
    params <- caribbean_10_model
    open <- !params@species_params$blocked_pred
    expect_true(any(open) && any(!open))
    expect_equal(getDiet(params, proportion = FALSE)[open, , ],
                 mizer_getDiet(params, proportion = FALSE)[open, , ],
                 ignore_attr = TRUE)
})

test_that("getDiet() equals mizer's method when no predator is blocked", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    sp <- species_params(params)
    sp$blocked_pred <- FALSE
    species_params(params) <- sp
    expect_equal(getDiet(params), mizer_getDiet(params), ignore_attr = TRUE)
})

test_that("getDiet() applies the refuge with default and positional arguments", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    explicit <- getDiet(params, proportion = TRUE, n = initialN(params),
                        n_pp = initialNResource(params),
                        n_other = initialNOther(params), t = 0)
    expect_equal(getDiet(params), explicit)
    expect_equal(getDiet(params, FALSE), getDiet(params, proportion = FALSE))
    expect_false(isTRUE(all.equal(getDiet(params), mizer_getDiet(params),
                                  ignore_attr = TRUE)))
})

test_that("getDiet() proportions sum to 1 wherever fish eat", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    total <- rowSums(getDiet(params), dims = 2)
    eats <- rowSums(getDiet(params, proportion = FALSE), dims = 2) > 0
    expect_equal(total[eats], rep(1, sum(eats)), ignore_attr = TRUE)
})

test_that("getDiet() keeps the diet of predator sizes that refuge hides completely", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@other_params$refuge_params$max_protect <- 1
    vul <- getVulnerable(params)
    expect_true(any(vul["predators", params@initial_n["predators", ] > 0] == 0))

    diet <- getDiet(params, proportion = FALSE)
    bal <- prey_balance(params, diet)
    expect_equal(bal["consumed", ], bal["lost", ], tolerance = 1e-8)
    expect_equal(rowSums(diet, dims = 2), expected_total(params),
                 tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("getDiet() uses the model's feeding level for satiating blocked predators", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    sp <- species_params(params)
    sp$satiation[sp$species == "predators"] <- TRUE
    species_params(params) <- sp
    expect_true(any(getFeedingLevel(params)["predators", ] > 0.01))

    diet <- getDiet(params, proportion = FALSE)
    bal <- prey_balance(params, diet)
    expect_equal(bal["consumed", ], bal["lost", ], tolerance = 1e-8)
    expect_equal(rowSums(diet, dims = 2), expected_total(params),
                 tolerance = 1e-8, ignore_attr = TRUE)
})

test_that("getDiet(sim) uses each saved time, including for refuge degradation", {
    data(caribbean_3_model, rubble_scale)
    params <- setDegradation(caribbean_3_model, deg_scale = rubble_scale,
                             bleach_time = 1, degrade = TRUE)
    sim <- project(params, t_max = 3, progress_bar = FALSE)
    diet_time <- getDiet(sim, proportion = FALSE)
    times <- as.numeric(dimnames(sim@n)$time)
    expect_equal(dim(diet_time)[1], length(times))
    expect_identical(names(dimnames(diet_time))[1], "time")

    for (i in seq_along(times)) {
        n <- array(sim@n[i, , ], dim = dim(sim@n)[2:3])
        dimnames(n) <- dimnames(sim@n)[2:3]
        n_other <- sim@n_other[i, ]
        names(n_other) <- dimnames(sim@n_other)$component
        expected <- getDiet(sim@params, proportion = FALSE, n = n,
                            n_pp = sim@n_pp[i, ], n_other = n_other,
                            t = times[i])
        expect_equal(diet_time[i, , , ], expected, ignore_attr = TRUE,
                     info = paste("time =", times[i]))
    }

    # After bleaching the refuge has changed, so the time matters
    last <- length(times)
    n <- array(sim@n[last, , ], dim = dim(sim@n)[2:3])
    dimnames(n) <- dimnames(sim@n)[2:3]
    n_other <- sim@n_other[last, ]
    names(n_other) <- dimnames(sim@n_other)$component
    at_t0 <- getDiet(sim@params, proportion = FALSE, n = n,
                     n_pp = sim@n_pp[last, ], n_other = n_other, t = 0)
    expect_false(isTRUE(all.equal(diet_time[last, , , ], at_t0,
                                  ignore_attr = TRUE)))
})

test_that("getDiet(sim) time_range and drop behave like mizer's", {
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 1, progress_bar = FALSE)
    kept <- getDiet(sim, time_range = 1)
    dropped <- getDiet(sim, time_range = 1, drop = TRUE)
    expect_identical(dim(kept)[1], 1L)
    expect_equal(dim(dropped), dim(kept)[-1])
})

test_that("plotDiet() shows the refuge-aware diet", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    plot_data <- plotDiet(params, species = "predators", return_data = TRUE)
    diet <- getDiet(params)
    row <- plot_data[plot_data$Prey == "predators", ][1, ]
    wi <- which.min(abs(params@w - row$w))
    expect_equal(row$Proportion, diet["predators", wi, "predators"],
                 tolerance = 1e-8, ignore_attr = TRUE)
    # mizer's own method shows more cannibalism, because it ignores refuge
    expect_lt(row$Proportion, mizer_getDiet(params)["predators", wi, "predators"])

    sim <- project(params, t_max = 1, progress_bar = FALSE)
    expect_s3_class(plotDiet(sim), "ggplot")
})

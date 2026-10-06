# mizer's own getDiet() methods, which ignore predation refuge
mizer_getDiet <- utils::getS3method("getDiet", "MizerParams")
mizer_getDiet_sim <- utils::getS3method("getDiet", "MizerSim")

# The model state at one saved time of a simulation
sim_state <- function(sim, i) {
    n <- array(sim@n[i, , ], dim = dim(sim@n)[2:3])
    dimnames(n) <- dimnames(sim@n)[2:3]
    n_other <- sim@n_other[i, ]
    names(n_other) <- dimnames(sim@n_other)$component
    list(n = n, n_pp = sim@n_pp[i, ], n_other = n_other,
         t = as.numeric(dimnames(sim@n)$time[[i]]))
}

initial_state <- function(params) {
    list(n = initialN(params), n_pp = initialNResource(params),
         n_other = initialNOther(params), t = 0)
}

# For each fish prey, the biomass all predators consume and the biomass that
# prey loses to predation. These must agree.
prey_balance <- function(params, diet, s = initial_state(params)) {
    pred_mort <- getPredMort(params, n = s$n, n_pp = s$n_pp,
                             n_other = s$n_other, time_range = s$t)
    dw <- rep(params@dw, each = nrow(s$n))
    sapply(params@species_params$species, function(prey) c(
        consumed = sum(diet[, , prey] * s$n * dw),
        lost = sum(pred_mort[prey, ] * s$n[prey, ] * params@w * params@dw)
    ))
}

# Total consumption equals encounter x (1 - feeding level), wherever fish exist
expected_total <- function(params, s = initial_state(params)) {
    # getEncounter() takes the time as `t`; the others take `time_range`
    encounter <- getEncounter(params, n = s$n, n_pp = s$n_pp,
                              n_other = s$n_other, t = s$t)
    feeding_level <- getFeedingLevel(params, n = s$n, n_pp = s$n_pp,
                                     n_other = s$n_other, time_range = s$t)
    total <- encounter * (1 - feeding_level)
    total[s$n <= 0] <- 0
    total
}

# Largest relative difference of any element. expect_equal() averages the
# difference over all elements, which can hide a single wrong size.
expect_each_equal <- function(actual, expected, info = NULL) {
    actual <- as.vector(actual)
    expected <- as.vector(expected)
    scale <- pmax(abs(expected), 1e-14 * max(abs(expected)))
    expect_lt(max(abs(actual - expected) / scale), 1e-8, label = info)
}

expect_diet_consistent <- function(params, diet, s = initial_state(params),
                                   info = NULL) {
    bal <- prey_balance(params, diet, s)
    expect_each_equal(bal["consumed", ], bal["lost", ], info = info)
    expect_each_equal(rowSums(diet, dims = 2), expected_total(params, s),
                      info = info)
}

test_that("getDiet() matches the model's consumption in both bundled models", {
    data(caribbean_3_model, caribbean_10_model)
    for (params in list(caribbean_3_model, caribbean_10_model)) {
        expect_diet_consistent(params, getDiet(params, proportion = FALSE))
    }
})

test_that("mizer's own getDiet() method overstates prey that use refuge", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    bal <- prey_balance(params, mizer_getDiet(params, proportion = FALSE))
    # predators use refuge, inverts don't
    expect_gt(bal["consumed", "predators"], 1.1 * bal["lost", "predators"])
    expect_equal(bal["consumed", "inverts"], bal["lost", "inverts"],
                 tolerance = 1e-8)
})

test_that("getDiet() matches the model at every time of a degrading simulation", {
    data(caribbean_3_model, rubble_scale)
    params <- setDegradation(caribbean_3_model, deg_scale = rubble_scale,
                             bleach_time = 1, degrade = TRUE)
    sim <- project(params, t_max = 3, progress_bar = FALSE)
    diet_time <- getDiet(sim, proportion = FALSE)
    times <- as.numeric(dimnames(sim@n)$time)
    expect_equal(dim(diet_time)[1], length(times))
    expect_identical(names(dimnames(diet_time))[1], "time")
    for (i in seq_along(times)) {
        expect_diet_consistent(sim@params, diet_time[i, , , ],
                               sim_state(sim, i),
                               info = paste("time =", times[i]))
    }

    # After bleaching the refuge has changed, so the time matters
    s <- sim_state(sim, length(times))
    at_t0 <- getDiet(sim@params, proportion = FALSE, n = s$n, n_pp = s$n_pp,
                     n_other = s$n_other, t = 0)
    expect_false(isTRUE(all.equal(diet_time[length(times), , , ], at_t0,
                                  ignore_attr = TRUE)))
})

test_that("getDiet() uses the given algae and detritus for the feeding level", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    s <- initial_state(params)
    s$n_other <- lapply(s$n_other, function(x) 3 * x)
    diet <- getDiet(params, proportion = FALSE, n_other = s$n_other)
    expect_diet_consistent(params, diet, s)
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

test_that("getDiet() applies the refuge with default, named and positional arguments", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    explicit <- getDiet(params, proportion = TRUE, n = initialN(params),
                        n_pp = initialNResource(params),
                        n_other = initialNOther(params), t = 0)
    expect_equal(getDiet(params), explicit)
    expect_equal(getDiet(params, FALSE), getDiet(params, proportion = FALSE))
    expect_equal(getDiet(params, FALSE, initialN(params), initialNResource(params)),
                 getDiet(params, proportion = FALSE))
    expect_false(isTRUE(all.equal(getDiet(params), mizer_getDiet(params),
                                  ignore_attr = TRUE)))
})

test_that("getDiet() runs an extension stacked above mizerReef once", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    # Defined here, S3 dispatch finds it without registering it globally
    getDiet.testStackedExt <- function(object, proportion = TRUE, ...) {
        2 * NextMethod()
    }
    stacked <- params
    class(stacked) <- c("testStackedExt", class(params))
    expect_equal(getDiet(stacked, proportion = FALSE),
                 2 * getDiet(params, proportion = FALSE))
    expect_equal(getDiet(stacked, proportion = FALSE, n = initialN(params),
                         n_pp = initialNResource(params),
                         n_other = initialNOther(params)),
                 2 * getDiet(params, proportion = FALSE))
})

test_that("getDiet() keeps t when an extension above passes positional arguments", {
    data(caribbean_3_model, rubble_scale)
    params <- setDegradation(caribbean_3_model, deg_scale = rubble_scale,
                             bleach_time = 1, degrade = TRUE)
    # The refuge changes at t = 1, so a wrong t would show
    expect_false(isTRUE(all.equal(getDiet(params), getDiet(params, t = 1))))
    # An extension above using the same named-NextMethod() pattern
    getDiet.testUpperExt <- function(object, proportion = TRUE,
                                     n = initialN(object),
                                     n_pp = initialNResource(object),
                                     n_other = initialNOther(object), ...) {
        NextMethod(proportion = proportion, n = n, n_pp = n_pp,
                   n_other = n_other)
    }
    upper <- params
    class(upper) <- c("testUpperExt", class(params))
    expect_equal(getDiet(upper, TRUE, initialN(params)), getDiet(params))
})

test_that("getDiet() asks for t rather than time_range on a params object", {
    data(caribbean_3_model)
    expect_error(getDiet(caribbean_3_model, time_range = 3), "single time as `t`")
})

test_that("getDiet() works with NaN abundances, as mizer's method does", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- initialN(params)
    n["herbivores", 50:60] <- NaN
    expect_no_error(mizer_getDiet(params, n = n))
    expect_no_error(getDiet(params, n = n))
    expect_no_error(getDiet(params, proportion = FALSE, n = n))
})

test_that("getDiet() equals mizer's method when nothing hides in refuge", {
    data(caribbean_3_model)
    params <- newRefuge(caribbean_3_model, new_method = "noncomplex")
    expect_true(all(getVulnerable(params) == 1))
    expect_equal(getDiet(params), mizer_getDiet(params), ignore_attr = TRUE)
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
    hidden <- params@initial_n["predators", ] > 0 & vul["predators", ] == 0
    expect_true(any(hidden))
    expect_diet_consistent(params, getDiet(params, proportion = FALSE))

    # Also for abundances so small that n * 1e-20 underflows to 0
    tiny <- which(hidden)[1]
    params@initial_n["predators", tiny] <- 1e-306
    diet <- getDiet(params, proportion = FALSE)
    expect_each_equal(rowSums(diet, dims = 2), expected_total(params))
})

test_that("getDiet() uses the model's feeding level for satiating blocked predators", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    sp <- species_params(params)
    sp$satiation[sp$species == "predators"] <- TRUE
    species_params(params) <- sp
    expect_true(any(getFeedingLevel(params)["predators", ] > 0.01))
    expect_diet_consistent(params, getDiet(params, proportion = FALSE))
})

test_that("getDiet(sim) selects times like mizer's method", {
    data(caribbean_3_model)
    sim <- project(caribbean_3_model, t_max = 2, t_save = 1, progress_bar = FALSE)
    # Without degradation the time doesn't change the diet, so mizer's
    # MizerSim method (which drops the time) gives the same values
    for (args in list(list(), list(time_range = 1), list(time_range = c(1, 2)),
                      list(time_range = 1, drop = TRUE),
                      list(proportion = FALSE, time_range = 2))) {
        expect_equal(do.call(getDiet, c(list(sim), args)),
                     do.call(mizer_getDiet_sim, c(list(sim), args)),
                     info = deparse(args))
    }
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

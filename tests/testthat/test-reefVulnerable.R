test_that("reefVulnerable equals 1 - refuge for the static methods (sigmoidal/binned/noncomplex)", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other

    for (method in list(
        list(m = "sigmoidal", mp = data.frame(L_refuge = 5, prop_protect = 0.6)),
        list(m = "binned", mp = data.frame(start_L = c(0, 10, 20), end_L = c(10, 20, 30), prop_protect = c(0.8, 0.5, 0.2))),
        list(m = "noncomplex", mp = NULL)
    )) {
        p <- getRefuge(setRefuge(params, method = method$m, method_params = method$mp))
        result <- reefVulnerable(p, n = n, n_pp = n_pp, n_other = n_other, t = 0)
        expected <- 1 - p@other_params$refuge_params$refuge
        expect_equal(unname(result), unname(expected), info = method$m)
    }
})

test_that("reefVulnerable computes the documented competitive-method formula with dummy fish bins", {
    # caribbean_3_model uses the competitive method. Hand-derive competitor
    # density and refuge proportions per bin directly from n/dw/bin.id/
    # refuge_density, independently of reefVulnerable()'s own code.
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other

    mp <- params@other_params$refuge_params$method_params
    bin_id <- params@other_params$refuge_params$bin.id
    max_protect <- params@other_params$refuge_params$max_protect
    tau <- params@other_params$refuge_params$tau
    refuge_user <- params@species_params$refuge_user
    no_w <- length(params@w)
    rd <- mp$refuge_density

    refuge <- matrix(0, nrow = 3, ncol = no_w)
    refuge_users <- params@species_params$species[refuge_user]
    for (k in seq_along(rd)) {
        idx <- bin_id[[k]]
        bin_fish <- sweep(n, 2, 1:no_w %in% idx, "*")
        cd <- sum(bin_fish[refuge_users, , drop = FALSE] %*% params@dw)
        refuge[, idx] <- if (cd == 0) max_protect else tau * rd[k] / cd
    }
    refuge[refuge > max_protect] <- max_protect
    expected <- 1 - (refuge_user * refuge)

    result <- reefVulnerable(params, n = n, n_pp = n_pp, n_other = n_other, t = 0)
    expect_equal(unname(result), unname(expected))
})

test_that("reefVulnerable computes the documented competitive-method formula with species-specific bins", {
    # caribbean_3_model's species share identical a/b, so it can't
    # distinguish species-specific from dummy bins; caribbean_10_species
    # has genuinely different a/b per species.
    data(caribbean_10_species)
    data(caribbean_10_interaction)
    sp <- caribbean_10_species
    base <- suppressMessages(newMultispeciesParams(
        species_params = sp, interaction = caribbean_10_interaction,
        min_w_pp = NA, w_pp_cutoff = 1, n = 0.75, p = 0.75
    ))
    base@other_params$refuge_params <- list()
    method_params <- data.frame(
        start_L = c(0, 10, 20), end_L = c(10, 20, 30), refuge_density = c(1, 0.5, 0.1)
    )
    params <- getRefuge(
        setRefuge(base, method = "competitive", method_params = method_params, use_dummy_fish_bins = FALSE),
        use_dummy_fish_bins = FALSE
    )
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other

    bin_id <- params@other_params$refuge_params$bin.id
    max_protect <- params@other_params$refuge_params$max_protect
    tau <- params@other_params$refuge_params$tau
    refuge_user <- params@species_params$refuge_user
    no_w <- length(params@w)
    no_sp <- nrow(sp)
    rd <- method_params$refuge_density
    refuge_user_idx <- which(refuge_user)

    refuge <- matrix(0, nrow = no_sp, ncol = no_w)
    for (i in seq_len(no_sp)) {
        for (k in seq_along(rd)) {
            idx <- bin_id[[paste0("sp", i, "_bin", k)]]
            bin_fish <- sweep(n, 2, 1:no_w %in% idx, "*")
            cd <- sum(bin_fish[refuge_user_idx, , drop = FALSE] %*% params@dw)
            refuge[i, idx] <- if (cd == 0) max_protect else tau * rd[k] / cd
        }
    }
    refuge[refuge > max_protect] <- max_protect
    expected <- 1 - (refuge_user * refuge)

    result <- reefVulnerable(params, n = n, n_pp = n_pp, n_other = n_other, t = 0)
    expect_equal(unname(result), unname(expected))
    # Confirm species genuinely differ (guards against the fixed bug where
    # species-specific bins silently collapsed to a shared profile).
    expect_false(isTRUE(all.equal(result[1, ], result[3, ])))
})

test_that("reefVulnerable gives max_protect within bin coverage when there are no competitors, and full vulnerability outside any bin", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    zero_n <- n
    zero_n[c("predators", "herbivores"), ] <- 0 # zero out all refuge-using species

    result <- reefVulnerable(params,
        n = zero_n, n_pp = params@initial_n_pp,
        n_other = params@initial_n_other, t = 0
    )

    bin_id <- params@other_params$refuge_params$bin.id
    covered_idx <- sort(unique(unlist(bin_id)))
    max_protect <- params@other_params$refuge_params$max_protect
    refuge_user <- params@species_params$refuge_user
    no_w <- length(params@w)

    expected <- matrix(1, nrow = 3, ncol = no_w) # uncovered sizes: refuge = 0, fully vulnerable
    for (i in 1:3) expected[i, covered_idx] <- 1 - refuge_user[i] * max_protect

    expect_equal(unname(result), unname(expected))
})

test_that("reefVulnerable with degrade = TRUE and new_rd = NULL matches passing reefDegrade()'s output explicitly", {
    data(caribbean_3_model)
    data(rubble_scale)
    params <- setDegradation(caribbean_3_model, deg_scale = rubble_scale, bleach_time = 2, degrade = TRUE)
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other

    auto <- reefVulnerable(params, n = n, n_pp = n_pp, n_other = n_other, t = 3)
    new_rd <- reefDegrade(params, n = n, n_pp = n_pp, n_other = n_other, t = 3)
    explicit <- reefVulnerable(params, n = n, n_pp = n_pp, n_other = n_other, t = 3, new_rd = new_rd)

    expect_equal(auto, explicit)
    # And it should differ from the pre-bleaching (t < t_bleach) profile.
    pre_bleach <- reefVulnerable(params, n = n, n_pp = n_pp, n_other = n_other, t = 0)
    expect_false(isTRUE(all.equal(auto, pre_bleach)))
})

test_that("reefVulnerable defaults a missing use_dummy_fish_bins to the dummy (TRUE) branch, matching getRefuge()", {
    # Regression test: reefVulnerable() used to default a missing/NULL
    # use_dummy_fish_bins to the species-specific branch (via isTRUE()),
    # opposite of getRefuge()'s isFALSE()-based default -- silently
    # producing an all-zero refuge matrix (100% vulnerability) whenever
    # bin.id was actually built the dummy way, which is exactly the state
    # of the bundled caribbean_3_model before this was fixed (it predates
    # the use_dummy_fish_bins parameter entirely).
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@other_params$refuge_params$use_dummy_fish_bins <- NULL
    expect_null(params@other_params$refuge_params$use_dummy_fish_bins)

    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    result <- reefVulnerable(params, n = n, n_pp = n_pp, n_other = n_other, t = 0)

    explicit_true <- params
    explicit_true@other_params$refuge_params$use_dummy_fish_bins <- TRUE
    expected <- reefVulnerable(explicit_true, n = n, n_pp = n_pp, n_other = n_other, t = 0)

    expect_equal(result, expected)
    expect_false(all(result == 1)) # the bug made every entry exactly 1
})

test_that("reefVulnerable warns when use_dummy_fish_bins = FALSE but bin.id doesn't have species-specific keys", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    # bin.id was built under the dummy convention (numeric keys); forcing
    # use_dummy_fish_bins = FALSE creates exactly the mismatch this guards.
    params@other_params$refuge_params$use_dummy_fish_bins <- FALSE

    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    expect_warning(
        result <- reefVulnerable(params, n = n, n_pp = n_pp, n_other = n_other, t = 0),
        "no predation refuge is being applied"
    )
    expect_true(all(result == 1))
})

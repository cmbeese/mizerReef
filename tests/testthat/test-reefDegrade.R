test_that("reefDegrade returns the current refuge density for the competitive method when degrade is FALSE", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_identical(params@other_params$refuge_params$method, "competitive")
    expect_false(isTRUE(params@other_params$refuge_params$degrade))

    result <- reefDegrade(
        params,
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other,
        t = 1
    )
    expect_identical(result, params@other_params$refuge_params$method_params$refuge_density)
})

test_that("reefDegrade returns NULL for a non-competitive method", {
    data(caribbean_3_model)
    data(tuning_profile)
    params <- newRefuge(caribbean_3_model,
        new_method = "binned",
        new_method_params = tuning_profile
    )
    result <- reefDegrade(
        params,
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other,
        t = 1
    )
    expect_null(result)
})

# Numeric coverage below hand-computes expected refuge densities
# independently of reefDegrade()'s own recursion, directly from the
# documented formula (deg_scale column 1 = bleaching year, columns 2+ =
# post-bleaching years 1, 2, 3...; at each step the column's scaling factors
# multiply the *previous* step's refuge density -- so, since the baseline is
# never mutated between calls, this is mathematically a cumulative product
# of deg_scale columns 1..k applied fresh to the baseline every time).
degrade_expected <- function(deg_scale, baseline_rd, t_bleach, t) {
    if (t < t_bleach) {
        return(baseline_rd)
    }
    rd <- baseline_rd
    deg_duration <- ncol(deg_scale)
    for (yr in t_bleach:t) {
        col <- if (yr == t_bleach) 1 else (yr - t_bleach + 1)
        if (col > deg_duration) next # exhausted: stays at last computed value
        rd <- deg_scale[, col] * rd
    }
    rd
}

test_that("reefDegrade matches hand-computed values at the bleach year, for all three built-in trajectories", {
    data(caribbean_3_model)
    data(rubble_scale)
    data(algae_scale)
    data(recovery_scale)
    params0 <- caribbean_3_model
    baseline_rd <- params0@other_params$refuge_params$method_params$refuge_density
    n <- params0@initial_n
    n_pp <- params0@initial_n_pp
    n_other <- params0@initial_n_other

    for (traj in list(rubble_scale, algae_scale, recovery_scale)) {
        params <- setDegradation(params0,
            deg_scale = traj, bleach_time = 2, degrade = TRUE
        )
        result <- reefDegrade(params, n = n, n_pp = n_pp, n_other = n_other, t = 2)
        expected <- degrade_expected(traj, baseline_rd, t_bleach = 2, t = 2)
        expect_equal(result, expected)
    }
})

test_that("reefDegrade matches hand-computed values across several post-bleaching years", {
    data(caribbean_3_model)
    data(rubble_scale)
    params0 <- caribbean_3_model
    baseline_rd <- params0@other_params$refuge_params$method_params$refuge_density
    n <- params0@initial_n
    n_pp <- params0@initial_n_pp
    n_other <- params0@initial_n_other

    params <- setDegradation(params0,
        deg_scale = rubble_scale, bleach_time = 2, degrade = TRUE
    )

    for (t in c(3, 4, 6, 10)) {
        result <- reefDegrade(params, n = n, n_pp = n_pp, n_other = n_other, t = t)
        expected <- degrade_expected(rubble_scale, baseline_rd, t_bleach = 2, t = t)
        expect_equal(result, expected, info = paste("t =", t))
    }
})

test_that("reefDegrade holds refuge density fixed once the deg_scale trajectory is exhausted", {
    data(caribbean_3_model)
    data(rubble_scale)
    params0 <- caribbean_3_model
    baseline_rd <- params0@other_params$refuge_params$method_params$refuge_density
    n <- params0@initial_n
    n_pp <- params0@initial_n_pp
    n_other <- params0@initial_n_other
    t_bleach <- 2
    deg_duration <- ncol(rubble_scale)

    params <- setDegradation(params0,
        deg_scale = rubble_scale, bleach_time = t_bleach, degrade = TRUE
    )

    t_last_col <- t_bleach + deg_duration - 1
    rd_last_col <- reefDegrade(params, n = n, n_pp = n_pp, n_other = n_other, t = t_last_col)
    rd_exhausted <- reefDegrade(params, n = n, n_pp = n_pp, n_other = n_other, t = t_last_col + 1)
    rd_further <- reefDegrade(params, n = n, n_pp = n_pp, n_other = n_other, t = t_last_col + 5)

    expected <- degrade_expected(rubble_scale, baseline_rd, t_bleach, t_last_col)
    expect_equal(rd_last_col, expected)
    expect_equal(rd_exhausted, rd_last_col)
    expect_equal(rd_further, rd_last_col)
})
test_that("reefPredMort computes the documented formula sum_i pred_rate_i * V_i * theta_i", {
    # Independent vectorized reimplementation of the documented formula
    # mu_p.i(w_p) = sum_j pred_rate_j(w_p) V_ji(w_p) theta_ji, using outer()
    # instead of reefPredMort()'s own row-replication loop.
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    pred_rate <- getPredRate(params)
    vulnerable <- reefVulnerable(params, n = n, n_pp = n_pp, n_other = n_other, t = 0)

    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    idx_sp <- (length(params@w_full) - no_w + 1):length(params@w_full)
    pr <- pred_rate[, idx_sp, drop = FALSE]
    int <- params@interaction
    blocked_pred <- params@species_params$blocked_pred == TRUE

    expected <- matrix(0, nrow = no_sp, ncol = no_w)
    for (i in 1:no_sp) {
        V_i <- if (blocked_pred[i]) vulnerable else matrix(1, no_sp, no_w)
        expected <- expected + V_i * outer(int[i, ], pr[i, ])
    }

    result <- reefPredMort(params, n, n_pp, n_other, t = 0, pred_rate = pred_rate)
    expect_equal(unname(result), unname(expected))
})

test_that("reefPredMort reduces to mizerPredMort() when no species is blocked_pred", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@species_params$blocked_pred <- FALSE
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    pred_rate <- getPredRate(params)

    result <- reefPredMort(params, n, n_pp, n_other, t = 0, pred_rate = pred_rate)
    expected <- mizer::mizerPredMort(params, n, n_pp, n_other, t = 0, pred_rate = pred_rate)
    expect_equal(unname(result), unname(expected))
})

test_that("reefPredMort respects interaction matrix orientation: zeroing a prey column zeroes that prey's mortality", {
    # Mirrors mizer's own "interaction is right way round" test pattern.
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@interaction[, "herbivores"] <- 0
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    pred_rate <- getPredRate(params)

    result <- reefPredMort(params, n, n_pp, n_other, t = 0, pred_rate = pred_rate)
    expect_true(all(result["herbivores", ] == 0))
})

test_that("projectPredMort.mizerReef dispatch matches reefPredMort called directly", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    pred_rate <- getPredRate(params)

    via_dispatch <- projectPredMort(params, n, n_pp, n_other, t = 0, pred_rate = pred_rate)
    via_direct <- reefPredMort(params, n, n_pp, n_other, t = 0, pred_rate = pred_rate)
    expect_equal(unname(via_dispatch), unname(via_direct))
})

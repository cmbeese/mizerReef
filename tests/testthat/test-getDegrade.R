test_that("getDegrade runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(getDegrade(params,
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other), NA)
})

test_that("getDegrade(params) matches reefDegrade() called directly at the same time", {
    data(caribbean_3_model)
    data(rubble_scale)
    params <- setDegradation(caribbean_3_model,
        deg_scale = rubble_scale, bleach_time = 2, degrade = TRUE
    )
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other

    for (t in c(0, 1, 2, 3, 5)) {
        via_getDegrade <- getDegrade(params,
            n = n, n_pp = n_pp, n_other = n_other, time_range = t
        )
        via_reefDegrade <- reefDegrade(params, n = n, n_pp = n_pp, n_other = n_other, t = t)
        expect_equal(via_getDegrade, via_reefDegrade, info = paste("t =", t))
    }
})

test_that("getDegrade(sim) matches reefDegrade() at each saved timestep", {
    data(caribbean_3_model)
    data(rubble_scale)
    params <- setDegradation(caribbean_3_model,
        deg_scale = rubble_scale, bleach_time = 2, degrade = TRUE
    )
    sim <- project(params, t_max = 5, progress_bar = FALSE)

    deg_time <- getDegrade(sim)
    times <- as.numeric(dimnames(sim@n)$time)
    expect_equal(dim(deg_time), c(length(times), 10))

    for (i in seq_along(times)) {
        t <- times[i]
        n <- array(sim@n[i, , ], dim = dim(sim@n)[2:3])
        dimnames(n) <- dimnames(sim@n)[2:3]
        n_other <- sim@n_other[i, ]
        names(n_other) <- dimnames(sim@n_other)$component
        expected <- reefDegrade(sim@params,
            n = n, n_pp = sim@n_pp[i, ], n_other = n_other, t = t
        )
        expect_equal(unname(deg_time[i, ]), unname(expected), info = paste("time =", t))
    }
})

test_that("getDegrade(sim) drop argument controls whether the time dimension is kept", {
    data(caribbean_3_model)
    data(rubble_scale)
    params <- setDegradation(caribbean_3_model,
        deg_scale = rubble_scale, bleach_time = 2, degrade = TRUE
    )
    # A single-timestep simulation gives a time dimension of length 1, so
    # drop = TRUE/FALSE produce different shapes -- a meaningful check that
    # the argument is honoured, not just that both calls run.
    sim <- project(params, t_max = 1, progress_bar = FALSE)

    dropped <- getDegrade(sim, time_range = 0, drop = TRUE)
    kept <- getDegrade(sim, time_range = 0, drop = FALSE)

    expect_null(dim(dropped))
    expect_identical(dim(kept), c(1L, 10L))
    expect_identical(names(dimnames(kept))[1], "time")
    expect_equal(unname(dropped), unname(kept[1, ]))
})

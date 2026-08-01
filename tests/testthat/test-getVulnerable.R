test_that("getVulnerable(params) matches reefVulnerable() called directly at the same time", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other

    via_getVulnerable <- getVulnerable(params, n = n, n_pp = n_pp, n_other = n_other)
    via_direct <- reefVulnerable(params, n = n, n_pp = n_pp, n_other = n_other, t = 0)
    expect_equal(unname(via_getVulnerable), unname(via_direct))
})

test_that("getVulnerable(params) sets dimnames to match params@metab", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- getVulnerable(params,
        n = params@initial_n, n_pp = params@initial_n_pp, n_other = params@initial_n_other
    )
    expect_equal(dimnames(result), dimnames(params@metab))
})

test_that("getVulnerable(sim) matches reefVulnerable() at each saved timestep", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    sim <- project(params, t_max = 3, progress_bar = FALSE)

    vul_time <- getVulnerable(sim)
    times <- as.numeric(dimnames(sim@n)$time)
    expect_equal(dim(vul_time), c(length(times), dim(sim@n)[2:3]))

    for (i in seq_along(times)) {
        t <- times[i]
        n <- array(sim@n[i, , ], dim = dim(sim@n)[2:3])
        dimnames(n) <- dimnames(sim@n)[2:3]
        n_other <- sim@n_other[i, ]
        names(n_other) <- dimnames(sim@n_other)$component
        expected <- reefVulnerable(sim@params, n = n, n_pp = sim@n_pp[i, ], n_other = n_other, t = t)
        expect_equal(unname(vul_time[i, , ]), unname(expected), info = paste("time =", t))
    }
})

test_that("getVulnerable(sim) drop argument controls whether the time dimension is kept", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    sim <- project(params, t_max = 1, progress_bar = FALSE)

    dropped <- getVulnerable(sim, time_range = 0, drop = TRUE)
    kept <- getVulnerable(sim, time_range = 0, drop = FALSE)

    expect_equal(dim(dropped), dim(sim@n)[2:3])
    expect_identical(dim(kept), c(1L, dim(sim@n)[2:3]))
    expect_identical(names(dimnames(kept))[1], "time")
    expect_equal(unname(dropped), unname(kept[1, , ]))
})

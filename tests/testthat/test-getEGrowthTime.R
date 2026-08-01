test_that("getEGrowthTime(params) matches directly calling the params' own EGrowth rate function", {
    # mizerReef leaves EGrowth at mizer's stock rate function (only the
    # five reef-specific rates -- Rates/Encounter/FeedingLevel/PredMort/
    # Mort -- are overridden); confirm getEGrowthTime() really does
    # dispatch to whatever @rates_funcs$EGrowth names, rather than
    # hardcoding or duplicating the formula.
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    f <- get(params@rates_funcs$EGrowth)

    expected <- f(params,
        n = n, n_pp = n_pp, n_other = n_other, t = 0,
        e_repro = getERepro(params, n = n, n_pp = n_pp, n_other = n_other, t = 0),
        e = getEReproAndGrowth(params, n = n, n_pp = n_pp, n_other = n_other, t = 0)
    )
    dimnames(expected) <- dimnames(params@metab)

    result <- getEGrowthTime(params, n = n, n_pp = n_pp, n_other = n_other)
    expect_equal(result, expected)
})

test_that("getEGrowthTime(params) has the same dimnames as params@metab", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- getEGrowthTime(params,
        n = params@initial_n, n_pp = params@initial_n_pp,
        n_other = params@initial_n_other
    )
    expect_equal(dimnames(result), dimnames(params@metab))
})

test_that("getEGrowthTime(sim) at a given time matches getEGrowthTime(params) at that time's abundances", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    sim <- project(params, t_max = 2, dt = 1, progress_bar = FALSE)

    result_sim <- getEGrowthTime(sim, time_range = 1)
    result_params <- getEGrowthTime(sim@params,
        n = sim@n["1", , ], n_pp = sim@n_pp["1", ],
        n_other = as.list(sim@n_other["1", ])
    )
    expect_equal(as.numeric(result_sim), as.numeric(result_params))
})

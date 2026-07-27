test_that("reefRates' individual rates match calling each underlying function directly with the same chained inputs", {
    # reefRates() orchestrates reefDegrade()/reefVulnerable()/reefEncounter()/
    # reefFeedingLevel()/reefPredMort()/reefMort() (plus the standard mizer
    # rate functions passed in via rates_fns) in sequence. Confirm each
    # returned rate matches an independent call to that same function, fed
    # with the previous rate in the chain -- an internal-consistency audit
    # rather than a re-derivation of reefRates()'s own code.
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    effort <- rep(0, nrow(params@species_params))
    rates_fns <- list(
        EReproAndGrowth = mizerEReproAndGrowth,
        ERepro = mizerERepro,
        EGrowth = mizerEGrowth,
        PredRate = mizerPredRate,
        FMort = mizerFMort,
        RDI = mizerRDI,
        RDD = BevertonHoltRDD,
        ResourceMort = mizerResourceMort
    )

    r <- reefRates(params,
        n = n, n_pp = n_pp, n_other = n_other, t = 0,
        effort = effort, rates_fns = rates_fns
    )

    expect_equal(r$degrade, reefDegrade(params, n, n_pp, n_other, t = 0))
    expect_equal(r$vulnerable, reefVulnerable(params, n, n_pp, n_other, t = 0, new_rd = r$degrade))
    expect_equal(r$encounter, reefEncounter(params, n, n_pp, n_other, t = 0, vulnerable = r$vulnerable))
    expect_equal(r$feeding_level, reefFeedingLevel(params, n, n_pp, n_other, t = 0, encounter = r$encounter))
    expect_equal(
        r$pred_rate,
        mizerPredRate(params, n, n_pp, n_other, feeding_level = r$feeding_level, vulnerable = r$vulnerable, t = 0)
    )
    expect_equal(r$pred_mort, reefPredMort(params, n, n_pp, n_other, t = 0, pred_rate = r$pred_rate))
    expect_equal(
        r$f_mort,
        mizerFMort(params, n, n_pp, n_other, effort = effort, t = 0, e_growth = r$e_growth, pred_mort = r$pred_mort)
    )
    expect_equal(r$mort, reefMort(params, n, n_pp, n_other, t = 0, f_mort = r$f_mort, pred_mort = r$pred_mort))
})

test_that("reefRates returns all documented rate components", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    rates_fns <- list(
        EReproAndGrowth = mizerEReproAndGrowth,
        ERepro = mizerERepro,
        EGrowth = mizerEGrowth,
        PredRate = mizerPredRate,
        FMort = mizerFMort,
        RDI = mizerRDI,
        RDD = BevertonHoltRDD,
        ResourceMort = mizerResourceMort
    )
    r <- reefRates(params,
        n = params@initial_n, n_pp = params@initial_n_pp, n_other = params@initial_n_other,
        t = 0, effort = rep(0, nrow(params@species_params)), rates_fns = rates_fns
    )
    expect_setequal(
        names(r),
        c(
            "degrade", "vulnerable", "encounter", "feeding_level", "e", "e_repro",
            "e_growth", "pred_rate", "pred_mort", "f_mort", "mort", "rdi", "rdd",
            "resource_mort"
        )
    )
})

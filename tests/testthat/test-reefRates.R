test_that("reefRates runs without error", {
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
    expect_error(reefRates(
        params,
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other,
        t = 0,
        effort = rep(0, nrow(params@species_params)),
        rates_fns = rates_fns
    ), NA)
})

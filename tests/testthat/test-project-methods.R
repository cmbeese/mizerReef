test_that("projectEncounter.mizerReef matches reefEncounter", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other

    via_dispatch <- projectEncounter(params, n, n_pp, n_other, t = 0)
    via_direct <- reefEncounter(params, n, n_pp, n_other, t = 0)
    # dimnames labelling can differ (mizer's base rate function names its
    # dimnames, reefEncounter() does not); only the values matter here
    expect_equal(unname(via_dispatch), unname(via_direct))
})

test_that("projectPredMort.mizerReef matches reefPredMort", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    pred_rate <- getPredRate(params)

    via_dispatch <- projectPredMort(params, n, n_pp, n_other, t = 0,
                                    pred_rate = pred_rate)
    via_direct <- reefPredMort(params, n, n_pp, n_other, t = 0,
                               pred_rate = pred_rate)
    expect_equal(unname(via_dispatch), unname(via_direct))
})

test_that("projectMort.mizerReef matches reefMort", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    f_mort <- getFMort(params)
    pred_mort <- getPredMort(params)

    via_dispatch <- projectMort(params, n, n_pp, n_other, t = 0,
                                f_mort = f_mort, pred_mort = pred_mort)
    via_direct <- reefMort(params, n, n_pp, n_other, t = 0,
                           f_mort = f_mort, pred_mort = pred_mort)
    expect_equal(via_dispatch, via_direct)
})

test_that("project() gives identical results when no species is blocked_pred", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    # Direct slot assignment, not species_params(params)$blocked_pred <- FALSE:
    # mizer 3.2.0's species_params<-.MizerParams only writes changed values
    # into the given_species_params tracking table, and since blocked_pred
    # was never tracked there for this object, species whose value didn't
    # change (already FALSE) get silently left as NA instead of FALSE.
    params@species_params$blocked_pred <- FALSE

    sim <- project(params, t_max = 2, progress_bar = FALSE)
    expect_false(anyNA(sim@n))
})

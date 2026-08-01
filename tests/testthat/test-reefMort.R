test_that("reefMort equals mizerMort() + reefSenMort() when include_sen_mort = TRUE", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_true(isTRUE(params@other_params$include_sen_mort))
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    f_mort <- getFMort(params)
    pred_mort <- getPredMort(params)

    expected <- mizer::mizerMort(params, n, n_pp, n_other, t = 0, f_mort = f_mort, pred_mort = pred_mort) +
        reefSenMort(params)
    result <- reefMort(params, n, n_pp, n_other, t = 0, f_mort = f_mort, pred_mort = pred_mort)
    expect_equal(unname(result), unname(expected))
})

test_that("reefMort equals mizerMort() alone when include_sen_mort = FALSE", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@other_params$include_sen_mort <- FALSE
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    f_mort <- getFMort(params)
    pred_mort <- getPredMort(params)

    expected <- mizer::mizerMort(params, n, n_pp, n_other, t = 0, f_mort = f_mort, pred_mort = pred_mort)
    result <- reefMort(params, n, n_pp, n_other, t = 0, f_mort = f_mort, pred_mort = pred_mort)
    expect_equal(unname(result), unname(expected))
})

test_that("projectMort.mizerReef dispatch matches reefMort called directly", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    f_mort <- getFMort(params)
    pred_mort <- getPredMort(params)

    via_dispatch <- projectMort(params, n, n_pp, n_other, t = 0, f_mort = f_mort, pred_mort = pred_mort)
    via_direct <- reefMort(params, n, n_pp, n_other, t = 0, f_mort = f_mort, pred_mort = pred_mort)
    expect_equal(via_dispatch, via_direct)
})

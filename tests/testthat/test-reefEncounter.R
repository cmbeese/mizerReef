test_that("reefEncounter reduces to mizerEncounter() when no species is blocked_pred", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@species_params$blocked_pred <- FALSE
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other

    result <- reefEncounter(params, n, n_pp, n_other, t = 0)
    expected <- mizer::mizerEncounter(params, n, n_pp, n_other, t = 0)
    expect_equal(unname(result), unname(expected))
})

test_that("reefEncounter's unblocked-predator rows match standard mizerEncounter() even when other species are blocked", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_true(any(params@species_params$blocked_pred))
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other

    result <- reefEncounter(params, n, n_pp, n_other, t = 0)
    standard <- mizer::mizerEncounter(params, n, n_pp, n_other, t = 0)
    good_pred <- which(params@species_params$blocked_pred == FALSE)
    expect_equal(unname(result[good_pred, ]), unname(standard[good_pred, ]))
})

test_that("reefEncounter's blocked-predator rows match mizerEncounter() computed with vulnerability-reduced prey abundance", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    vulnerable <- reefVulnerable(params, n = n, n_pp = n_pp, n_other = n_other, t = 0)

    result <- reefEncounter(params, n, n_pp, n_other, t = 0)
    standard_with_vul_n <- mizer::mizerEncounter(params, n = vulnerable * n, n_pp = n_pp, n_other = n_other, t = 0)
    blocked_pred <- which(params@species_params$blocked_pred == TRUE)
    expect_equal(unname(result[blocked_pred, ]), unname(standard_with_vul_n[blocked_pred, ]))
})

test_that("projectEncounter.mizerReef dispatch matches reefEncounter called directly", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other

    via_dispatch <- projectEncounter(params, n, n_pp, n_other, t = 0)
    via_direct <- reefEncounter(params, n, n_pp, n_other, t = 0)
    expect_equal(unname(via_dispatch), unname(via_direct))
})

test_that("reefSteady returns a valid mizerReef object with finite, non-negative abundances", {
    # Mirrors mizer's own test-steady.R pattern (all(is.finite(...)),
    # all(... >= 0)) applied to the mizerReef-extended object.
    data(caribbean_3_model)
    result <- reefSteady(caribbean_3_model, progress_bar = FALSE)
    expect_s4_class(result, "mizerReef")
    expect_true(all(is.finite(result@initial_n)))
    expect_true(all(result@initial_n >= 0))
})

test_that("reefSteady reaches a genuine fixed point: one more t_per step barely moves the state", {
    # Stronger than trusting reefSteady()'s own internal convergence check:
    # project one further t_per (1.5 years) step from the result and confirm
    # mizer's own distanceSSLogN() (used by projectToSteady() to *define*
    # convergence) stays below the same default tolerance (tol = 0.1*dt).
    data(caribbean_3_model)
    result <- reefSteady(caribbean_3_model, progress_bar = FALSE)

    previous <- list(
        n = result@initial_n, n_pp = result@initial_n_pp,
        n_other = result@initial_n_other
    )
    sim <- project(result, t_max = 1.5, t_save = 1.5, dt = 0.1, progress_bar = FALSE)
    last <- dim(sim@n)[1]
    current <- list(
        n = array(sim@n[last, , ], dim = dim(sim@n)[2:3], dimnames = dimnames(sim@n)[2:3]),
        n_pp = sim@n_pp[last, ],
        n_other = list(
            algae = sim@n_other[last, "algae"],
            detritus = sim@n_other[last, "detritus"]
        )
    )
    d <- distanceSSLogN(result, current, previous)
    expect_lt(d, 0.1 * 0.1)
})

test_that("reefSteady preserves reproduction_level by default", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    old_level <- mizer::getReproductionLevel(params)
    result <- reefSteady(params, progress_bar = FALSE)
    expect_equal(unname(mizer::getReproductionLevel(result)), unname(old_level))
})

test_that("reefSteady preserves R_max when preserve = 'R_max'", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    old_R_max <- params@species_params$R_max
    result <- reefSteady(params, preserve = "R_max", progress_bar = FALSE)
    expect_equal(unname(result@species_params$R_max), unname(old_R_max))
})

test_that("reefSteady preserves erepro when preserve = 'erepro', except where mizer must clamp it upward", {
    # setBevertonHolt() can bump erepro above the requested value (with a
    # warning) when the requested value can't sustain the given reproduction
    # rate -- the same benign warning documented for matchReadyGrowth() in
    # test-matchReefGrowth.R. So the meaningful assertion is "unchanged or
    # increased", not exact equality.
    data(caribbean_3_model)
    params <- caribbean_3_model
    old_erepro <- params@species_params$erepro
    result <- suppressWarnings(reefSteady(params, preserve = "erepro", progress_bar = FALSE))
    expect_true(all(result@species_params$erepro >= old_erepro - 1e-10))
})

test_that("reefSteady retunes algae/detritus to a genuine steady state when new_refuge = FALSE", {
    data(caribbean_3_model)
    result <- reefSteady(caribbean_3_model, progress_bar = FALSE)

    P_A <- sum(getAlgaeProduction(result))
    c_A <- algae_consumption(result, n = result@initial_n, rates = getRates(result))
    expect_equal(P_A - c_A * algae_biomass(result), 0)

    P_D <- sum(getDetritusProduction(result))
    c_D <- detritus_consumption(result, n = result@initial_n, rates = getRates(result))
    expect_equal(P_D - c_D * detritus_biomass(result), 0)
})

test_that("reefSteady does not retune algae/detritus when new_refuge = TRUE", {
    data(caribbean_3_model)
    mp <- data.frame(L_refuge = 0.5, prop_protect = 0.3)
    new_refuge_params <- newRefuge(caribbean_3_model,
        new_refuge = TRUE, new_method = "sigmoidal", new_method_params = mp
    )
    old_growth <- new_refuge_params@other_params$algae$growth
    old_external <- new_refuge_params@other_params$detritus$external

    result <- reefSteady(new_refuge_params, progress_bar = FALSE)
    expect_equal(result@other_params$algae$growth, old_growth)
    expect_equal(result@other_params$detritus$external, old_external)
})

test_that("reefSteady dispatches to tuneUR_cc when use_UR_cc is TRUE", {
    data(caribbean_3_model)
    cc_params <- setURcapacity(caribbean_3_model, cap = 1.5)
    result <- suppressWarnings(reefSteady(cc_params, progress_bar = FALSE))

    ba <- algae_biomass(result)
    ka <- result@other_params$algae$capacity
    c_A <- algae_consumption(result, n = result@initial_n, rates = getRates(result))
    P_A <- sum(getAlgaeProduction(result))
    expect_equal(P_A * (1 - ba / ka) - c_A * ba, 0)
})

test_that("reefSteady with return_sim = TRUE returns a MizerSim-like object wrapping the tuned params", {
    data(caribbean_3_model)
    old_growth <- caribbean_3_model@other_params$algae$growth
    sim <- reefSteady(caribbean_3_model, return_sim = TRUE, progress_bar = FALSE)
    expect_s4_class(sim@params, "mizerReef")

    # algae_growth is a fixed input, left unchanged by tuning -- only the
    # algae biomass is retuned to the resulting steady state.
    expect_equal(sim@params@other_params$algae$growth, old_growth)
    P_A <- sum(getAlgaeProduction(sim@params))
    c_A <- algae_consumption(sim@params, n = sim@params@initial_n, rates = getRates(sim@params))
    expect_equal(P_A - c_A * algae_biomass(sim@params), 0)
})

test_that("loading mizerReef leaves mizer::steady() working for non-reef models", {
    # mizerReef used to install reefSteady() over mizer::steady() with
    # assignInNamespace(), which replaced mizer's S3 generic outright: every
    # non-reef model in the session then hit reefSteady()'s reef-only code
    # (`params@other_params$new_refuge` is NULL there, so `if (NULL == FALSE)`
    # errored), and no other extension's steady() method could dispatch.
    expect_true(grepl("UseMethod", paste(deparse(body(mizer::steady)),
                                         collapse = " ")))
    # NS_params is nowhere near steady after 3 years, so setBevertonHolt()
    # warns about erepro; that is beside the point being tested here.
    reef_free <- suppressWarnings(suppressMessages(
        mizer::steady(mizer::NS_params, t_max = 3, progress_bar = FALSE)
    ))
    expect_s4_class(reef_free, "MizerParams")
    expect_false(is(reef_free, "mizerReef"))
})

test_that("steady() and tuneSteadyState() dispatch to reefSteady() for reef models", {
    data(caribbean_3_model)
    expected <- reefSteady(caribbean_3_model, progress_bar = FALSE)

    via_steady <- mizer::steady(caribbean_3_model, progress_bar = FALSE)
    via_tune <- mizer::tuneSteadyState(caribbean_3_model, progress_bar = FALSE)

    expect_s4_class(via_steady, "mizerReef")
    expect_s4_class(via_tune, "mizerReef")
    expect_equal(via_steady@initial_n, expected@initial_n)
    expect_equal(via_tune@initial_n, expected@initial_n)
    expect_equal(via_steady@initial_n_other, expected@initial_n_other)
    expect_equal(via_tune@initial_n_other, expected@initial_n_other)
})

test_that("tuneSteadyState() on a reef model refuses the Newton solver", {
    data(caribbean_3_model)
    expect_error(mizer::tuneSteadyState(caribbean_3_model, solver = "newton"),
                 "supports only")
})

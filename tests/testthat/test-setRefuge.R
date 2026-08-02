test_that("setRefuge stores method, method_params and default scalar parameters", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    method_params <- data.frame(
        start_L = 1,
        end_L = 10,
        prop_protect = 0.5
    )
    result <- setRefuge(params, method = "binned", method_params = method_params)

    expect_identical(result@other_params$refuge_params$method, "binned")
    expect_equal(result@other_params$refuge_params$method_params, method_params)
    # Documented defaults
    expect_equal(result@other_params$refuge_params$a_bar, 0.025)
    expect_equal(result@other_params$refuge_params$b_bar, 3)
    expect_equal(result@other_params$refuge_params$w_settle, 0.1)
    expect_equal(result@other_params$refuge_params$max_protect, 0.98)
    expect_equal(result@other_params$refuge_params$tau, 1)
    expect_true(result@other_params$refuge_params$use_dummy_fish_bins)
})

test_that("setRefuge converts list method_params to a data frame identically to a data frame input", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    mp_df <- data.frame(start_L = 1, end_L = 10, prop_protect = 0.5)
    mp_list <- list(start_L = 1, end_L = 10, prop_protect = 0.5)

    from_df <- setRefuge(params, method = "binned", method_params = mp_df)
    from_list <- setRefuge(params, method = "binned", method_params = mp_list)

    expect_equal(
        from_df@other_params$refuge_params$method_params,
        from_list@other_params$refuge_params$method_params
    )
})

test_that("setRefuge stores user-supplied scalar parameters instead of defaults", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    method_params <- data.frame(start_L = 1, end_L = 10, prop_protect = 0.5)

    result <- setRefuge(params,
        method = "binned", method_params = method_params,
        a_bar = 0.03, b_bar = 3.1, w_settle = 0.2, max_protect = 0.9, tau = 0.5
    )

    expect_equal(result@other_params$refuge_params$a_bar, 0.03)
    expect_equal(result@other_params$refuge_params$b_bar, 3.1)
    expect_equal(result@other_params$refuge_params$w_settle, 0.2)
    expect_equal(result@other_params$refuge_params$max_protect, 0.9)
    expect_equal(result@other_params$refuge_params$tau, 0.5)
})

test_that("setRefuge defaults slope to 100 for the sigmoidal method when not supplied", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- setRefuge(params,
        method = "sigmoidal",
        method_params = data.frame(L_refuge = 0.5, prop_protect = 0.3)
    )
    expect_equal(result@other_params$refuge_params$method_params$slope, 100)
})

test_that("setRefuge fills missing refuge_user/blocked_pred with FALSE and warns", {
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    sp <- caribbean_3_species
    params <- suppressMessages(newMultispeciesParams(
        species_params = sp, interaction = caribbean_3_interaction,
        min_w_pp = NA, w_pp_cutoff = 1, n = 0.75, p = 0.75
    ))
    params@other_params$refuge_params <- list()

    sp_no_refuge_user <- params
    sp_no_refuge_user@species_params$refuge_user <- NULL
    expect_warning(
        result <- setRefuge(sp_no_refuge_user, method = "noncomplex"),
        "no species use refuge"
    )
    expect_true(all(!result@species_params$refuge_user))

    sp_no_blocked_pred <- params
    sp_no_blocked_pred@species_params$blocked_pred <- NULL
    expect_warning(
        result2 <- setRefuge(sp_no_blocked_pred, method = "noncomplex"),
        "all predators can access prey within refuge"
    )
    expect_true(all(!result2@species_params$blocked_pred))
})

test_that("setRefuge computes the documented default satiation rule and warns", {
    # satiation default = !eats_other_species & eats_detritus & !eats_algae,
    # i.e. TRUE only for pure detritivores -- derived from the interaction
    # matrix and the interaction_algae/interaction_detritus columns --
    # independently hand-computed here from caribbean_3_interaction/species,
    # not by re-deriving from setRefuge's own code.
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    sp <- caribbean_3_species
    params <- suppressMessages(newMultispeciesParams(
        species_params = sp, interaction = caribbean_3_interaction,
        min_w_pp = NA, w_pp_cutoff = 1, n = 0.75, p = 0.75
    ))
    params@other_params$refuge_params <- list()
    params@species_params$satiation <- NULL

    eats_other_species <- rowSums(caribbean_3_interaction) > 0
    eats_algae <- sp[["interaction_algae"]] > 0
    eats_detritus <- sp[["interaction_detritus"]] > 0
    expected_satiation <- !eats_other_species & eats_detritus & !eats_algae

    expect_warning(
        result <- setRefuge(params, method = "noncomplex"),
        "default is TRUE only for detritivores"
    )
    expect_equal(unname(result@species_params$satiation), unname(expected_satiation))
    # For caribbean_3_species specifically, the tightened default reproduces
    # the bundled caribbean_3_model's own explicit, thesis-matched values
    # exactly (predators/herbivores = FALSE, inverts = TRUE) -- confirming
    # the heuristic no longer conflates algae-grazing herbivores with
    # detritivores the way the old "any resource consumer" rule did.
    expect_equal(
        unname(result@species_params$satiation),
        c(FALSE, FALSE, TRUE)
    )
})

test_that("setRefuge fills missing a/b with a_bar/b_bar and warns", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@species_params$a <- NULL
    params@species_params$b <- NULL
    method_params <- data.frame(start_L = 1, end_L = 10, prop_protect = 0.5)

    expect_warning(
        result <- setRefuge(params,
            method = "binned", method_params = method_params,
            a_bar = 0.04, b_bar = 3.2
        ),
        "set to average values"
    )
    expect_true(all(result@species_params$a == 0.04))
    expect_true(all(result@species_params$b == 3.2))
})

test_that("setRefuge errors on an invalid method, missing method_params, and bad values", {
    data(caribbean_3_model)
    params <- caribbean_3_model

    expect_error(setRefuge(params, method = "bogus"), "Method must be")
    expect_error(setRefuge(params, method = "binned"), "must provide method specific parameters")
    expect_error(
        setRefuge(params,
            method = "binned",
            method_params = data.frame(start_L = 10, end_L = 1, prop_protect = 0.5)
        ),
        "start lengths must be less than bin end lengths"
    )
    expect_error(
        setRefuge(params,
            method = "binned",
            method_params = data.frame(start_L = 1, end_L = 10, prop_protect = 1.5)
        ),
        "prop_protect should be a proportion between 0 and 1"
    )
    expect_error(
        setRefuge(params,
            method = "binned",
            method_params = data.frame(start_L = 1, end_L = 10, prop_protect = 0.5),
            a_bar = -1
        ),
        "a_bar must be non-negative"
    )
})

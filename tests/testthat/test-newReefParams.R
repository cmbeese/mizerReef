test_that("newReefParams returns a mizerReef object for the binned method", {
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    data(tuning_profile)
    result <- suppressMessages(newReefParams(
        species_params = caribbean_3_species,
        interaction = caribbean_3_interaction,
        method = "binned",
        method_params = tuning_profile
    ))
    expect_s4_class(result, "mizerReef")
    expect_identical(result@other_params$refuge_params$method, "binned")
})

test_that("newReefParams returns a mizerReef object for the competitive method", {
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    data(caribbean_3_model)
    competitive_mp <- caribbean_3_model@other_params$refuge_params$method_params
    result <- suppressMessages(newReefParams(
        species_params = caribbean_3_species,
        interaction = caribbean_3_interaction,
        method = "competitive",
        method_params = competitive_mp
    ))
    expect_s4_class(result, "mizerReef")
    expect_identical(result@other_params$refuge_params$method, "competitive")
})

test_that("newReefParams's refuge matrix matches setRefuge + getRefuge called directly on the base params", {
    # newReefParams() internally calls newMultispeciesParams() then
    # setRefuge()/getRefuge() -- reproduce that low-level path by hand and
    # confirm the high-level wrapper doesn't diverge from it.
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    data(tuning_profile)
    result <- suppressMessages(newReefParams(
        species_params = caribbean_3_species,
        interaction = caribbean_3_interaction,
        method = "binned",
        method_params = tuning_profile
    ))

    base <- suppressMessages(newMultispeciesParams(
        species_params = caribbean_3_species,
        interaction = caribbean_3_interaction,
        min_w_pp = NA, w_pp_cutoff = 1, n = 0.75, p = 0.75
    ))
    base@other_params$refuge_params <- list()
    base@other_params$algae <- list()
    base@other_params$detritus <- list()
    direct <- getRefuge(setRefuge(base, method = "binned", method_params = tuning_profile))

    expect_equal(
        result@other_params$refuge_params$refuge,
        direct@other_params$refuge_params$refuge,
        ignore_attr = TRUE
    )
    expect_equal(
        result@other_params$refuge_params$refuge_lengths,
        direct@other_params$refuge_params$refuge_lengths,
        ignore_attr = TRUE
    )
})

test_that("newReefParams computes rho_algae/rho_detritus using the documented formula", {
    # rho = pmax(0, f0 * h / (1 - f0) - E) * interaction_{algae,detritus},
    # hand-computed independently from newReefParams()'s own code, using
    # mizer::getEncounter() on the plain (pre-reef) params object.
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    data(tuning_profile)
    sp <- caribbean_3_species

    base <- suppressMessages(newMultispeciesParams(
        species_params = sp, interaction = caribbean_3_interaction,
        min_w_pp = NA, w_pp_cutoff = 1, n = 0.75, p = 0.75
    ))
    f0 <- mizer::set_species_param_default(base@species_params, "f0", 0.7)$f0
    E <- mizer::getEncounter(base)[, length(base@w)] / (base@w[length(base@w)]^0.75)
    expected_rho_algae <- pmax(0, f0 * base@species_params$h / (1 - f0) - E) *
        sp$interaction_algae
    expected_rho_detritus <- pmax(0, f0 * base@species_params$h / (1 - f0) - E) *
        sp$interaction_detritus

    result <- suppressMessages(newReefParams(
        species_params = sp,
        interaction = caribbean_3_interaction,
        method = "binned",
        method_params = tuning_profile
    ))

    expect_equal(unname(result@species_params$rho_algae), unname(expected_rho_algae))
    expect_equal(unname(result@species_params$rho_detritus), unname(expected_rho_detritus))
})

test_that("newReefParams computes external mortality using the allometric z0pre formula when include_ext_mort = FALSE", {
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    data(tuning_profile)
    result <- suppressMessages(newReefParams(
        species_params = caribbean_3_species,
        interaction = caribbean_3_interaction,
        method = "binned",
        method_params = tuning_profile,
        include_ext_mort = FALSE,
        z0pre = 0.15,
        n = 0.75
    ))
    expected_mort <- outer(rep(0.15, 3), result@w^0.25)
    expect_equal(mizer::ext_mort(result), expected_mort, ignore_attr = TRUE)
})

test_that("newReefParams computes external mortality using nat_mort = 0.2 default when include_ext_mort = TRUE", {
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    data(tuning_profile)
    result <- suppressMessages(newReefParams(
        species_params = caribbean_3_species,
        interaction = caribbean_3_interaction,
        method = "binned",
        method_params = tuning_profile,
        include_ext_mort = TRUE,
        n = 0.75
    ))
    expected_mort <- outer(rep(0.2, 3), result@w^0.25)
    expect_equal(mizer::ext_mort(result), expected_mort, ignore_attr = TRUE)
})

test_that("newReefParams stores degradation parameters when degrade = TRUE", {
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    data(tuning_profile)
    data(rubble_scale)
    result <- suppressMessages(newReefParams(
        species_params = caribbean_3_species,
        interaction = caribbean_3_interaction,
        method = "binned",
        method_params = tuning_profile,
        degrade = TRUE,
        deg_scale = rubble_scale,
        bleach_time = 3
    ))
    expect_true(result@other_params$refuge_params$degrade)
    expect_equal(result@other_params$refuge_params$t_bleach, 3)
    expect_equal(result@other_params$refuge_params$deg_scale, as.matrix(rubble_scale), ignore_attr = TRUE)
})

test_that("newReefParams sets degrade = FALSE by default", {
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    data(tuning_profile)
    result <- suppressMessages(newReefParams(
        species_params = caribbean_3_species,
        interaction = caribbean_3_interaction,
        method = "binned",
        method_params = tuning_profile
    ))
    expect_false(result@other_params$refuge_params$degrade)
})

test_that("newReefParams errors when method is missing", {
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    data(tuning_profile)
    expect_error(
        newReefParams(
            species_params = caribbean_3_species,
            interaction = caribbean_3_interaction,
            method_params = tuning_profile
        ),
        "method"
    )
})

test_that("newReefParams populates other_params$algae/detritus consistently with mizer::getComponent(), for both use_UR_cc settings", {
    # Regression check: other_params$algae/other_params$detritus is the
    # single, mizer-canonical storage location for these components
    # (mizer::setComponent()'s/getComponent()'s contract), populated once at
    # construction via setComponent()'s component_params argument. This
    # checks it is actually populated and that mizer::getComponent() -- a
    # real, exported mizer function any user or extension might call --
    # sees the same values newReefParams() itself just set, for both
    # use_UR_cc branches.
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    data(tuning_profile)

    for (use_UR_cc in c(FALSE, TRUE)) {
        result <- suppressMessages(newReefParams(
            species_params = caribbean_3_species,
            interaction = caribbean_3_interaction,
            method = "binned",
            method_params = tuning_profile,
            use_UR_cc = use_UR_cc
        ))

        expect_false(is.null(result@other_params$algae$growth))
        expect_equal(
            mizer::getComponent(result, "algae")$component_params$growth,
            result@other_params$algae$growth
        )
        expect_false(is.null(result@other_params$detritus$external))
        expect_equal(
            mizer::getComponent(result, "detritus")$component_params$external,
            result@other_params$detritus$external
        )
    }
})

test_that("newReefParams stores the algae/detritus consumption matrix under the correct field name", {
    # Regression check for the same migration-gap bug class as Finding 2:
    # newReefParams() used to write the consumption matrix to
    # other_params$algae_params$rho_algae / detritus_params$rho_detritus,
    # while every consumer (algae_components.R, detritus_components.R,
    # reef-components.R, reef-helpers.R) reads/writes the bare `rho` field
    # of the single, mizer-canonical other_params$algae/detritus location.
    # This only worked by accident of R's `$` partial name-matching (`rho`
    # uniquely prefix-matches `rho_algae`), which would silently break under
    # `[[` access or if a second field starting with `rho` were ever added.
    # Guard the exact field name on a freshly built object.
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    data(tuning_profile)

    result <- suppressMessages(newReefParams(
        species_params = caribbean_3_species,
        interaction = caribbean_3_interaction,
        method = "binned",
        method_params = tuning_profile
    ))

    a_names <- names(result@other_params$algae)
    d_names <- names(result@other_params$detritus)
    expect_true("rho" %in% a_names)
    expect_false("rho_algae" %in% a_names)
    expect_true("rho" %in% d_names)
    expect_false("rho_detritus" %in% d_names)
    expect_null(result@other_params$algae_params)
    expect_null(result@other_params$detritus_params)
})

test_that("caribbean_3_model's algae and detritus are genuinely at steady state", {
    # Regression check for Finding 2: caribbean_3_model's bundled
    # detritus_params external flux used to be stored under the stale field
    # name d_external (a leftover from an old S4-slot-then-revert migration),
    # so getDetritusProduction() -- which reads detritus_params$external
    # specifically -- silently dropped it. The bundled object's own on-disk
    # dB_D/dt was actually 120.1861, not zero, despite existing tuneUR tests
    # only ever checking a *fresh* tuneUR() call's output rather than the
    # shipped object's own state. This test checks the shipped object
    # directly, the same way test-caribbean_10_model.R does.
    data(caribbean_3_model)
    m <- with_calibration_senescence(caribbean_3_model)

    # Algae has a carrying capacity, detritus does not
    expect_equal(m@other_dynamics$algae, "algae_dynamics_cc")
    expect_equal(m@other_params$detritus$capacity, Inf)

    P_A <- sum(getAlgaeProduction(m))
    K_A <- m@other_params$algae$capacity
    B_A <- algae_biomass(m)
    c_A <- algae_consumption(m, n = m@initial_n, rates = getRates(m))
    expect_equal(P_A * (1 - B_A / K_A) - c_A * B_A, 0, tolerance = 1e-8)

    P_D <- sum(getDetritusProduction(m))
    c_D <- detritus_consumption(m, n = m@initial_n, rates = getRates(m))
    expect_equal(P_D - c_D * detritus_biomass(m), 0, tolerance = 1e-8)
})

test_that("caribbean_3_model's detritus external flux is stored under the correct field name", {
    # Direct guard against Finding 2's exact failure mode: the field must be
    # named `external` (what getDetritusProduction()/tuneUR()/tuneUR_cc()
    # actually read), not the stale `d_external` name.
    data(caribbean_3_model)
    dp_names <- names(caribbean_3_model@other_params$detritus)
    expect_true("external" %in% dp_names)
    expect_false("d_external" %in% dp_names)
    expect_false(is.null(caribbean_3_model@other_params$detritus$external))
})

test_that("caribbean_3_model's algae/detritus consumption matrix is stored under the correct field name", {
    # Regression check for the same migration-gap bug class as Finding 2:
    # newReefParams() used to write these matrices to other_params$
    # algae_params$rho_algae / detritus_params$rho_detritus, but every
    # consumer (algae_components.R, detritus_components.R,
    # reef-components.R, reef-helpers.R) reads/writes the bare `rho` field.
    # This only worked by accident of R's `$` partial name-matching (`rho`
    # uniquely prefix-matches `rho_algae`) - fragile, and would silently
    # break under `[[` access or if a second field starting with `rho` were
    # ever added. Guard the exact field name directly, bypassing `$`.
    data(caribbean_3_model)
    ap_names <- names(caribbean_3_model@other_params$algae)
    dp_names <- names(caribbean_3_model@other_params$detritus)
    expect_true("rho" %in% ap_names)
    expect_false("rho_algae" %in% ap_names)
    expect_true("rho" %in% dp_names)
    expect_false("rho_detritus" %in% dp_names)
})

test_that("caribbean_3_model no longer carries a duplicate algae_params/detritus_params structure", {
    # Regression check for the fix consolidating other_params$algae_params/
    # detritus_params onto the mizer-canonical other_params$algae/detritus:
    # only one structure should exist per resource.
    data(caribbean_3_model)
    expect_null(caribbean_3_model@other_params$algae_params)
    expect_null(caribbean_3_model@other_params$detritus_params)
})

test_that("caribbean_3_species and caribbean_3_interaction rebuild caribbean_3_model", {
    # The same call as inst/scripts/caribbean-3-calibration-clean.r, which
    # builds the bundled model from these tables. Only the initial
    # abundances differ: the script then projects to steady state.
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    data(karpata_refuge)
    data(caribbean_3_model)
    p <- newReefParams(
        species_params = caribbean_3_species,
        interaction = caribbean_3_interaction,
        method = "competitive", method_params = karpata_refuge,
        w_settle = 0.025 * 2.5^3, w_pp_cutoff = 10^-1.5, kappa = 0.1757,
        use_UR_cc = TRUE, algae_capacity = 902.8, detritus_capacity = Inf,
        external = 7.757, feces_decomp = 1, info_level = 0)
    m <- caribbean_3_model
    # The bundled model predates feces_decomp, which it treats as 1
    m@other_params$detritus$feces_decomp <- 1
    for (slot in c("search_vol", "intake_max", "metab", "mu_b", "psi",
                   "interaction", "rr_pp", "cc_pp")) {
        expect_equal(p[[slot]], m[[slot]], info = slot)
    }
    expect_equal(p@species_params, m@species_params)
    expect_equal(p@other_params$algae, m@other_params$algae)
    expect_mapequal(p@other_params$detritus, m@other_params$detritus)
    expect_equal(p@other_dynamics, m@other_dynamics)
})

test_that("caribbean_3_model is at steady state", {
    data(caribbean_3_model)
    sim <- project(with_calibration_senescence(caribbean_3_model),
                   t_max = 5, progress_bar = FALSE)
    b <- getBiomass(sim)
    expect_lt(max(abs(b[nrow(b), ] / b[1, ] - 1)), 1e-6)
})

test_that("caribbean_3_model matches its calibration targets", {
    # The targets in R/data.R and the vignette: observed biomasses above the
    # cutoff, Rogers et al. (2018)'s adult predator diet, predators maturing
    # at age 2, algae at 350 g/m^2 and a 3.85 h detritus lifetime.
    data(caribbean_3_model)
    m <- caribbean_3_model
    sp <- m@species_params
    expect_equal(getBiomass(m, use_cutoff = TRUE), sp$biomass_observed,
                 tolerance = 0.01, ignore_attr = TRUE)

    diet <- getDiet(m, proportion = FALSE)
    adult <- m@w >= sp$w_mat[1]
    eaten <- colSums(diet["predators", adult, ] * m@initial_n[1, adult] * m@dw[adult])
    share <- eaten[c("predators", "herbivores", "inverts")] / sum(eaten)
    expect_equal(unname(share), c(0.33, 0.17, 0.5), tolerance = 0.02)

    expect_equal(age_mat(m)[["predators"]], 2, tolerance = 0.01)
    expect_equal(algae_biomass(m), 350, tolerance = 0.01)
    expect_equal(detritus_lifetime(m) * 8760, 3.85, tolerance = 0.01)
})

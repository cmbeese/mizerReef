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
    m <- caribbean_3_model

    P_A <- sum(getAlgaeProduction(m))
    c_A <- algae_consumption(m, n = m@initial_n, rates = getRates(m))
    expect_equal(P_A - c_A * algae_biomass(m), 0, tolerance = 1e-8)

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

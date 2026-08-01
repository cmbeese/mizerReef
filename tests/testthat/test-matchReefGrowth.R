test_that("matchReefGrowth rescales algae/detritus rho by age_mat(params) / observed age_mat, for species with known age_mat", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    sp <- params@species_params
    sel <- !is.na(sp$age_mat)
    expect_true(any(sel)) # sanity check on the fixture

    factor <- mizer::age_mat(params)[sel] / sp$age_mat[sel]
    old_algae_rho <- params@other_params$algae$rho
    old_detritus_rho <- params@other_params$detritus$rho
    old_rho_algae <- params@species_params$rho_algae
    old_rho_detritus <- params@species_params$rho_detritus

    result <- matchReefGrowth(params)

    expect_equal(result@other_params$algae$rho[sel, ], old_algae_rho[sel, ] * factor)
    expect_equal(result@other_params$detritus$rho[sel, ], old_detritus_rho[sel, ] * factor)
    expect_equal(result@species_params$rho_algae[sel], old_rho_algae[sel] * factor)
    expect_equal(result@species_params$rho_detritus[sel], old_rho_detritus[sel] * factor)
})

test_that("matchReefGrowth leaves rho untouched for species with no known age_mat", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    sel <- is.na(params@species_params$age_mat)
    skip_if(!any(sel), "no species with unknown age_mat in caribbean_3_model")
    old_algae_rho <- params@other_params$algae$rho

    result <- matchReefGrowth(params)
    expect_equal(result@other_params$algae$rho[sel, ], old_algae_rho[sel, ])
})

test_that("matchReefGrowth restricted to a species subset only rescales that subset", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    known_idx <- which(!is.na(params@species_params$age_mat))
    target_idx <- known_idx[1]
    other_known_idx <- setdiff(known_idx, target_idx)
    skip_if(length(other_known_idx) == 0, "need at least two species with known age_mat")
    target <- params@species_params$species[target_idx]
    old_rho <- params@other_params$algae$rho

    result <- matchReefGrowth(params, species = target)
    expect_equal(result@other_params$algae$rho[other_known_idx, ], old_rho[other_known_idx, ])
})

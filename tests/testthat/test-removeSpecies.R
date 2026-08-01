test_that("removeSpecies drops the removed species' row from algae and detritus rho", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    sp_to_remove <- params@species_params$species[1]
    keep_idx <- which(params@species_params$species != sp_to_remove)
    old_algae_rho <- params@other_params$algae$rho
    old_detritus_rho <- params@other_params$detritus$rho

    result <- removeSpecies(params, sp_to_remove)

    expect_equal(nrow(result@other_params$algae$rho), nrow(params@species_params) - 1)
    expect_equal(result@other_params$algae$rho, old_algae_rho[keep_idx, , drop = FALSE])
    expect_equal(result@other_params$detritus$rho, old_detritus_rho[keep_idx, , drop = FALSE])
})

test_that("removeSpecies leaves the kept species' species_params rows unchanged", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    sp_to_remove <- params@species_params$species[1]
    kept <- setdiff(params@species_params$species, sp_to_remove)

    result <- removeSpecies(params, sp_to_remove)
    expect_equal(result@species_params$species, kept)
    expect_false(sp_to_remove %in% result@species_params$species)
})

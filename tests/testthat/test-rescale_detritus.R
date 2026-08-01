test_that("rescale_detritus multiplies biomass and divides rho/rho_detritus by the same factor", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    old_biomass <- detritus_biomass(params)
    old_rho <- params@other_params$detritus$rho
    old_rho_detritus <- params@species_params$rho_detritus

    result <- rescale_detritus(params, factor = 4)

    expect_equal(detritus_biomass(result), old_biomass * 4)
    expect_equal(result@other_params$detritus$rho, old_rho / 4)
    expect_equal(result@species_params$rho_detritus, old_rho_detritus / 4)
})

test_that("rescale_detritus keeps total detritus consumption (rate * biomass) unchanged", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    old_total <- detritus_consumption(params) * detritus_biomass(params)

    result <- rescale_detritus(params, factor = 3.7)
    new_total <- detritus_consumption(result) * detritus_biomass(result)

    expect_equal(new_total, old_total)
})

test_that("rescale_detritus does not touch decomposition proportions", {
    # rescale_detritus()'s own title says "without changing anything else"
    # -- confirm sen_decomp/ext_decomp (which govern getDetritusProduction()'s
    # decomp term, unrelated to biomass or rho) are genuinely untouched.
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- rescale_detritus(params, factor = 4)

    expect_equal(result@other_params$detritus$sen_decomp,
                 params@other_params$detritus$sen_decomp)
    expect_equal(result@other_params$detritus$ext_decomp,
                 params@other_params$detritus$ext_decomp)
})

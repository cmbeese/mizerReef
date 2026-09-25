test_that("scaleDownPlankton scales the plankton down and the search volume up", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- scaleDownPlankton(params, factor = 4)

    expect_equal(result@initial_n_pp, params@initial_n_pp / 4)
    expect_equal(result@cc_pp, params@cc_pp / 4)
    expect_equal(result@resource_params$kappa, params@resource_params$kappa / 4)
    expect_equal(result@search_vol, params@search_vol * 4)
    expect_equal(species_params(result)$gamma, species_params(params)$gamma * 4)
})

test_that("scaleDownPlankton leaves fish, algae and detritus unchanged", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- scaleDownPlankton(params, factor = 4)

    expect_equal(result@initial_n, params@initial_n)
    expect_equal(result@initial_n_other, params@initial_n_other)
    expect_equal(result@other_params$algae, params@other_params$algae)
    expect_equal(result@other_params$detritus, params@other_params$detritus)
    expect_equal(species_params(result)$rho_algae, species_params(params)$rho_algae)
    expect_equal(species_params(result)$rho_detritus, species_params(params)$rho_detritus)
})

test_that("scaleDownPlankton keeps each species' reproduction level", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- scaleDownPlankton(params, factor = 4)

    expect_equal(getReproductionLevel(result), getReproductionLevel(params))
})

test_that("scaleDownPlankton multiplies only the encounter rate with fish prey by the factor", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- scaleDownPlankton(params, factor = 4)

    # With no fish, what is left is the encounter with plankton, algae and
    # detritus and the external encounter, all of which should be unchanged.
    # The rate arrays carry their params object as an attribute, which does
    # differ, so compare only the rates themselves.
    no_fish <- initialN(params) * 0
    non_fish_before <- getEncounter(params, n = no_fish)
    non_fish_after <- getEncounter(result, n = no_fish)
    expect_equal(non_fish_after, non_fish_before, ignore_attr = "params")

    fish_before <- getEncounter(params) - non_fish_before
    fish_after <- getEncounter(result) - non_fish_after
    expect_equal(fish_after, fish_before * 4, ignore_attr = "params")
})

test_that("scaleDownPlankton leaves a recruitment function other than Beverton-Holt alone", {
    data(caribbean_3_model)
    params <- setReproduction(caribbean_3_model, RDD = "noRDD")
    params@species_params$R_max <- NULL
    result <- scaleDownPlankton(params, factor = 4)

    expect_equal(result@rates_funcs$RDD, "noRDD")
    expect_equal(result@species_params$erepro, params@species_params$erepro)
    expect_equal(result@search_vol, params@search_vol * 4)
})

test_that("scaleDownPlankton scales background species down with the plankton", {
    data(caribbean_3_model)
    params <- markBackground(caribbean_3_model, "inverts")
    result <- scaleDownPlankton(params, factor = 4)

    expect_equal(result@initial_n["inverts", ], params@initial_n["inverts", ] / 4)
    expect_equal(result@initial_n[c("predators", "herbivores"), ],
                 params@initial_n[c("predators", "herbivores"), ])
})

test_that("scaleDownPlankton with factor = 1 changes no rates", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- scaleDownPlankton(params, factor = 1)

    expect_equal(getEncounter(result), getEncounter(params), ignore_attr = "params")
    expect_equal(getEGrowth(result), getEGrowth(params), ignore_attr = "params")
    expect_equal(getMort(result), getMort(params), ignore_attr = "params")
})

test_that("scaleDownPlankton rejects non-positive factors", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(scaleDownPlankton(params, factor = 0))
    expect_error(scaleDownPlankton(params, factor = -2))
})

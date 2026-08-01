test_that("reefFeedingLevel matches the documented formula encounter / (encounter + intake_max) with satiation-adjusted intake_max", {
    # Following mizer's own test pattern for mizerFeedingLevel(): reconstruct
    # the formula independently rather than re-deriving reefFeedingLevel()'s
    # own code.
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    encounter <- reefEncounter(params, n, n_pp, n_other, t = 0)

    intake_max <- params@intake_max
    intake_max[params@species_params$satiation == FALSE] <- Inf
    expected <- encounter / (encounter + intake_max)
    expected[is.na(expected)] <- 0

    result <- reefFeedingLevel(params, n, n_pp, n_other, t = 0, encounter = encounter)
    expect_equal(unname(result), unname(expected))
})

test_that("reefFeedingLevel is exactly 0 for species with satiation = FALSE (unlimited intake)", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_true(any(params@species_params$satiation == FALSE))
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    encounter <- reefEncounter(params, n, n_pp, n_other, t = 0)

    result <- reefFeedingLevel(params, n, n_pp, n_other, t = 0, encounter = encounter)
    no_satiation <- params@species_params$satiation == FALSE
    expect_true(all(result[no_satiation, ] == 0))
})

test_that("reefFeedingLevel matches mizerFeedingLevel() directly for species with satiation = TRUE", {
    # For satiation = TRUE species, intake_max is untouched, so reefFeedingLevel()
    # should reduce to the standard mizer calculation for those rows.
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    encounter <- reefEncounter(params, n, n_pp, n_other, t = 0)

    result <- reefFeedingLevel(params, n, n_pp, n_other, t = 0, encounter = encounter)
    expected_mizer <- mizer::mizerFeedingLevel(params, n, n_pp, n_other, t = 0, encounter = encounter)
    satiation <- params@species_params$satiation == TRUE
    expect_equal(unname(result[satiation, ]), unname(expected_mizer[satiation, ]))
})

test_that("projectFeedingLevel.mizerReef dispatch matches reefFeedingLevel called directly", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    encounter <- reefEncounter(params, n, n_pp, n_other, t = 0)

    via_dispatch <- projectFeedingLevel(params, n, n_pp, n_other, t = 0, encounter = encounter)
    via_direct <- reefFeedingLevel(params, n, n_pp, n_other, t = 0, encounter = encounter)
    expect_equal(via_dispatch, via_direct)
})

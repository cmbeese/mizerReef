test_that("scaleReefAbundance multiplies all foreground species' abundance by a single factor", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    is_foreground <- !is.na(params@A)
    old_n <- params@initial_n

    result <- scaleReefAbundance(params, factor = 2)
    expect_equal(result@initial_n[is_foreground, ], old_n[is_foreground, ] * 2)
    # Background species (if any) are untouched.
    if (any(!is_foreground)) {
        expect_equal(result@initial_n[!is_foreground, ], old_n[!is_foreground, ])
    }
})

test_that("scaleReefAbundance applies a named factor only to the named species", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    sp <- params@species_params$species[1]
    old_n <- params@initial_n
    factor <- setNames(5, sp)

    result <- scaleReefAbundance(params, factor = factor)
    expect_equal(result@initial_n[sp, ], old_n[sp, ] * 5)

    other_sp <- setdiff(params@species_params$species, sp)
    expect_equal(result@initial_n[other_sp, ], old_n[other_sp, ])
})

test_that("scaleReefAbundance errors on an unknown species name", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(
        scaleReefAbundance(params, factor = setNames(2, "not_a_real_species")),
        "not_a_real_species"
    )
})

test_that("scaleReefAbundance rejects non-positive factors", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(scaleReefAbundance(params, factor = -1))
    expect_error(scaleReefAbundance(params, factor = 0))
})

test_that("reefSenMort computes the documented senescence formula", {
    # sen_mort[i,] = sen_prop * (log10(w)/log10(w_max[i]))^sen_curve, with
    # the ratio (not sen_prop * ratio) clipped to 0 before exponentiation --
    # hand-computed independently from params@w/w_max/ext_mort_params, not
    # by re-deriving reefSenMort()'s own code. This is the formula from the
    # original thesis chapter (Chapter3.tex eq-extmort): it approaches
    # exactly sen_prop as w -> w_max, which the clip-then-multiply order
    # used before this fix did not.
    data(caribbean_3_model)
    params <- caribbean_3_model
    mort_params <- params@other_params$ext_mort_params
    sen_prop <- mort_params$sen_prop
    sen_curve <- mort_params$sen_curve
    w <- params@w
    w_max <- params@species_params$w_max

    expected <- matrix(0, nrow = 3, ncol = length(w))
    for (i in 1:3) {
        ratio <- log10(w) / log10(w_max[i])
        ratio[ratio < 0] <- 0
        expected[i, ] <- sen_prop * ratio^sen_curve
    }

    result <- reefSenMort(params)
    expect_equal(unname(result), unname(expected))
})

test_that("reefSenMort approaches exactly sen_prop at w_max", {
    # Regression guard for the specific bug this fixed: sen_prop used to be
    # multiplied in before exponentiating, so the rate at w_max was
    # sen_prop^sen_curve (0.501 under the package defaults), not sen_prop
    # (0.1) as the thesis and roxygen document.
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- reefSenMort(params)
    sen_prop <- params@other_params$ext_mort_params$sen_prop

    for (i in 1:3) {
        expect_equal(result[i, ncol(result)], sen_prop)
    }
})

test_that("reefSenMort is zero at the smallest sizes and increases towards w_max", {
    # A structural sanity check independent of the exact formula: senescence
    # mortality should vanish for small individuals (log(w)/log(w_max) small,
    # clipped at 0) and be strictly increasing up to w_max.
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- reefSenMort(params)

    for (i in 1:3) {
        expect_equal(result[i, 1], 0)
        expect_true(all(diff(result[i, ]) >= 0))
    }
})

test_that("reefSenMort responds to setExtMortParams sen_prop/sen_curve overrides", {
    data(caribbean_3_model)
    params <- setExtMortParams(caribbean_3_model,
        ext_mort_params = data.frame(nat_mort = 0.2, sen_prop = 0.5, sen_curve = 2)
    )
    w <- params@w
    w_max <- params@species_params$w_max
    expected <- matrix(0, nrow = 3, ncol = length(w))
    for (i in 1:3) {
        ratio <- log10(w) / log10(w_max[i])
        ratio[ratio < 0] <- 0
        expected[i, ] <- 0.5 * ratio^2
    }

    result <- reefSenMort(params)
    expect_equal(unname(result), unname(expected))
})

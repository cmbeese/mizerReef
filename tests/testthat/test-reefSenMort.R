test_that("reefSenMort computes the documented senescence formula", {
    # sen_mort = sen_prop * (w / w_sen)^sen_curve for every species, the
    # senescence mortality of Rogers et al. (2014, 2018). Computed by hand
    # from params@w and ext_mort_params rather than by re-deriving
    # reefSenMort()'s own code.
    data(caribbean_3_model)
    params <- setExtMortParams(caribbean_3_model)
    mort_params <- params@other_params$ext_mort_params
    w <- params@w
    expected <- t(replicate(
        3, mort_params$sen_prop * (w / mort_params$w_sen)^mort_params$sen_curve))

    result <- reefSenMort(params)
    expect_equal(unname(result), unname(expected))
})

test_that("reefSenMort equals sen_prop at w_sen and increases with size", {
    # Rogers et al. (2018) Table S2: k_s = 0.1, m_s = 10^3 g, p_s = 0.3. The
    # size grid has no point at exactly 1000 g, so w_sen is put on the grid
    # point nearest 1 kg, where the rate must be exactly sen_prop.
    data(caribbean_3_model)
    w <- caribbean_3_model@w
    w_sen <- w[which.min(abs(log10(w) - 3))]
    params <- setExtMortParams(caribbean_3_model,
        ext_mort_params = list(nat_mort = 0.2, sen_prop = 0.1,
                               sen_curve = 0.3, w_sen = w_sen))
    result <- reefSenMort(params)
    expect_equal(unname(result[, w == w_sen]), rep(0.1, 3))
    # Increasing with size, and the same for every species
    for (i in 1:3) {
        expect_true(all(diff(result[i, ]) > 0))
        expect_equal(result[i, ], result[1, ])
    }
})

test_that("reefSenMort responds to setExtMortParams overrides", {
    data(caribbean_3_model)
    params <- setExtMortParams(caribbean_3_model,
        ext_mort_params = data.frame(nat_mort = 0.2, sen_prop = 0.5,
                                     sen_curve = 2, w_sen = 100)
    )
    w <- params@w
    result <- reefSenMort(params)
    expect_equal(unname(result[2, ]), 0.5 * (w / 100)^2)
})

test_that("reefSenMort uses w_sen = 1000 g for models without w_sen", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@other_params$ext_mort_params$w_sen <- NULL
    mort_params <- params@other_params$ext_mort_params
    expect_equal(unname(reefSenMort(params)[1, ]),
                 mort_params$sen_prop *
                     (params@w / 1000)^mort_params$sen_curve)
})

test_that("reefSenMort uses the power form for models without sen_form", {
    data(caribbean_3_model)
    params <- setExtMortParams(caribbean_3_model)
    expected <- reefSenMort(params)
    params@other_params$ext_mort_params$sen_form <- NULL
    expect_equal(reefSenMort(params), expected)
})

test_that("reefSenMort computes the deprecated log-ratio form", {
    # sen_prop * max(0, log10(w) / log10(w_max[i]))^sen_curve, the form
    # mizerReef used before 2.1.0, computed by hand per species
    data(caribbean_3_model)
    withr::local_options(lifecycle_verbosity = "quiet")
    params <- setExtMortParams(caribbean_3_model, sen_form = "log_ratio")
    mort_params <- params@other_params$ext_mort_params
    w <- params@w
    w_max <- params@species_params$w_max

    expected <- matrix(0, nrow = 3, ncol = length(w))
    for (i in 1:3) {
        ratio <- log10(w) / log10(w_max[i])
        ratio[ratio < 0] <- 0
        expected[i, ] <- mort_params$sen_prop * ratio^mort_params$sen_curve
    }
    result <- reefSenMort(params)
    expect_equal(unname(result), expected)

    # Zero below 1 g and sen_prop at w_max, whatever w_sen is
    for (i in 1:3) {
        expect_equal(result[i, 1], 0)
        expect_equal(result[i, ncol(result)], mort_params$sen_prop)
    }
    params@other_params$ext_mort_params$w_sen <- 1
    expect_equal(reefSenMort(params), result)
})

test_that("reefSenMort rejects an unknown stored sen_form", {
    data(caribbean_3_model)
    params <- setExtMortParams(caribbean_3_model)
    params@other_params$ext_mort_params$sen_form <- "log-ratio"
    expect_error(reefSenMort(params), "should be \"power\" or \"log_ratio\"")
})

test_that("reefSenMort returns a species x size array with dimnames", {
    data(caribbean_3_model)
    withr::local_options(lifecycle_verbosity = "quiet")
    for (form in c("power", "log_ratio")) {
        params <- setExtMortParams(caribbean_3_model, sen_form = form)
        expect_identical(dimnames(reefSenMort(params)),
                         dimnames(params@initial_n))
    }
})

test_that("reefSenMort errors when the model has no senescence settings", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@other_params$ext_mort_params <- NULL
    expect_error(reefSenMort(params), "no senescence settings")
})

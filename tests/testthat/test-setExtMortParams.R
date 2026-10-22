test_that("setExtMortParams defaults match its documented values", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- setExtMortParams(params)

    expect_equal(result@other_params$ext_mort_params,
                 list(nat_mort = 0.2, sen_prop = 0.1, sen_curve = 0.3,
                      w_sen = 1000, sen_form = "power"))
})

test_that("setExtMortParams stores a supplied named list, adding w_sen and sen_form", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    custom <- list(nat_mort = 0.5, sen_prop = 0.2, sen_curve = 0.6)

    result <- setExtMortParams(params, ext_mort_params = custom)
    expect_equal(result@other_params$ext_mort_params,
                 c(custom, list(w_sen = 1000, sen_form = "power")))

    custom$w_sen <- 500
    result <- setExtMortParams(params, ext_mort_params = custom)
    expect_equal(result@other_params$ext_mort_params,
                 c(custom, list(sen_form = "power")))
})

test_that("setExtMortParams normalises a data.frame/matrix input to the same list as a named list", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    as_list <- list(nat_mort = 0.5, sen_prop = 0.2, sen_curve = 0.6)
    as_df <- data.frame(nat_mort = 0.5, sen_prop = 0.2, sen_curve = 0.6)
    as_mat <- as.matrix(as_df)

    result_list <- setExtMortParams(params, ext_mort_params = as_list)
    result_df <- setExtMortParams(params, ext_mort_params = as_df)
    result_mat <- setExtMortParams(params, ext_mort_params = as_mat)

    expect_equal(result_list@other_params$ext_mort_params,
                 result_df@other_params$ext_mort_params)
    expect_equal(result_list@other_params$ext_mort_params,
                 result_mat@other_params$ext_mort_params)
})

test_that("setExtMortParams rejects negative values, missing required names and a bad w_sen", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(
        setExtMortParams(params, ext_mort_params = list(nat_mort = -1, sen_prop = 0.1, sen_curve = 0.3)),
        "nonnegative"
    )
    expect_error(
        setExtMortParams(params, ext_mort_params = list(sen_prop = 0.1, sen_curve = 0.3)),
        "nat_mort"
    )
    for (bad in list(0, NA_real_, Inf)) {
        expect_error(
            setExtMortParams(params, ext_mort_params = list(
                nat_mort = 0.2, sen_prop = 0.1, sen_curve = 0.3, w_sen = bad)),
            "positive number"
        )
    }
})

test_that("setExtMortParams rejects per-species values", {
    # The rates are the same for every species; a vector would be recycled
    # along the size axis by reefSenMort()
    data(caribbean_3_model)
    expect_error(
        setExtMortParams(caribbean_3_model, data.frame(
            nat_mort = 0.2, sen_prop = c(0.1, 0.2, 0.3), sen_curve = 0.3)),
        "single number"
    )
    expect_error(
        setExtMortParams(caribbean_3_model, list(
            nat_mort = 0.2, sen_prop = 0.1, sen_curve = 0.3,
            w_sen = c(500, 1000))),
        "single number"
    )
})

test_that("setExtMortParams accepts a named numeric vector", {
    data(caribbean_3_model)
    result <- setExtMortParams(caribbean_3_model,
        c(nat_mort = 0.5, sen_prop = 0.2, sen_curve = 0.6))
    expect_equal(result@other_params$ext_mort_params,
                 list(nat_mort = 0.5, sen_prop = 0.2, sen_curve = 0.6,
                      w_sen = 1000, sen_form = "power"))
    expect_error(setExtMortParams(caribbean_3_model, c(0.5, 0.2, 0.6)),
                 "named list")
})

test_that("setExtMortParams warns about and drops unknown names", {
    # A misspelt w_sen must not be taken as w_sen (partial matching) or
    # leave an empty entry behind
    data(caribbean_3_model)
    expect_warning(
        result <- setExtMortParams(caribbean_3_model, ext_mort_params = list(
            nat_mort = 0.2, sen_prop = 0.1, sen_curve = 0.3,
            w_senescence = 500)),
        "unknown external mortality parameters: w_senescence"
    )
    expect_equal(result@other_params$ext_mort_params,
                 list(nat_mort = 0.2, sen_prop = 0.1, sen_curve = 0.3,
                      w_sen = 1000, sen_form = "power"))
})

test_that("setExtMortParams rejects an unknown sen_form", {
    data(caribbean_3_model)
    expect_error(setExtMortParams(caribbean_3_model, sen_form = "linear"),
                 "should be one of")
})

test_that("sen_form = 'log_ratio' is deprecated but stored", {
    data(caribbean_3_model)
    lifecycle::expect_deprecated(
        result <- setExtMortParams(caribbean_3_model, sen_form = "log_ratio"),
        "log-ratio senescence form"
    )
    expect_equal(result@other_params$ext_mort_params$sen_form, "log_ratio")
})

test_that("the stored list can be passed back, keeping its sen_form unless sen_form is given", {
    data(caribbean_3_model)
    withr::local_options(lifecycle_verbosity = "quiet")
    old <- setExtMortParams(caribbean_3_model, sen_form = "log_ratio")
    stored <- old@other_params$ext_mort_params

    again <- setExtMortParams(old, ext_mort_params = stored)
    expect_equal(again@other_params$ext_mort_params, stored)

    switched <- setExtMortParams(old, ext_mort_params = stored,
                                 sen_form = "power")
    expect_equal(switched@other_params$ext_mort_params$sen_form, "power")
})

test_that("newReefParams passes sen_form to setExtMortParams", {
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    data(tuning_profile)
    build <- function(...) {
        suppressMessages(newReefParams(
            species_params = caribbean_3_species,
            interaction = caribbean_3_interaction,
            method = "binned", method_params = tuning_profile,
            info_level = 0, ...))
    }
    expect_equal(build()@other_params$ext_mort_params$sen_form, "power")
    lifecycle::expect_deprecated(p <- build(sen_form = "log_ratio"))
    expect_equal(p@other_params$ext_mort_params$sen_form, "log_ratio")

    # A stored list passed back keeps its form, as in setExtMortParams(),
    # unless sen_form is given
    stored <- p@other_params$ext_mort_params
    lifecycle::expect_deprecated(p2 <- build(ext_mort_params = stored))
    expect_equal(p2@other_params$ext_mort_params, stored)
    p3 <- build(ext_mort_params = stored, sen_form = "power")
    expect_equal(p3@other_params$ext_mort_params$sen_form, "power")
})

test_that("setExtMortParams defaults match its documented values", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- setExtMortParams(params)

    expect_equal(result@other_params$ext_mort_params,
                 list(nat_mort = 0.2, sen_prop = 0.1, sen_curve = 0.3))
})

test_that("setExtMortParams stores a supplied named list as-is", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    custom <- list(nat_mort = 0.5, sen_prop = 0.2, sen_curve = 0.6)

    result <- setExtMortParams(params, ext_mort_params = custom)
    expect_equal(result@other_params$ext_mort_params, custom)
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

test_that("setExtMortParams rejects negative values and missing required names", {
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
})

test_that("setDetritusParams defaults match its documented values", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- setDetritusParams(params)

    expect_equal(result@other_params$detritus$capacity, 1)
    expect_equal(result@other_params$detritus$sen_decomp, 0.8)
    expect_equal(result@other_params$detritus$ext_decomp, 0.2)
    expect_equal(result@other_params$detritus$external, 1)
    expect_false(result@other_params$use_UR_cc)
})

test_that("setDetritusParams stores the supplied values", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- setDetritusParams(params,
        detritus_capacity = 5, sen_decomp = 0.3,
        ext_decomp = 0.4, external = 99
    )

    expect_equal(result@other_params$detritus$capacity, 5)
    expect_equal(result@other_params$detritus$sen_decomp, 0.3)
    expect_equal(result@other_params$detritus$ext_decomp, 0.4)
    expect_equal(result@other_params$detritus$external, 99)
})

test_that("setDetritusParams rejects sen_decomp/ext_decomp outside [0, 1]", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(setDetritusParams(params, sen_decomp = 1.5), "proportion")
    expect_error(setDetritusParams(params, ext_decomp = -0.1), "proportion")
})

test_that("setDetritusParams warns and zeroes interaction_detritus when the column is missing", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@species_params$interaction_detritus <- NULL

    expect_warning(result <- setDetritusParams(params), "interaction_detritus")
    expect_true(all(result@species_params$interaction_detritus == 0))
})

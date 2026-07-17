test_that("reefDegrade returns the current refuge density for the competitive method when degrade is FALSE", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_identical(params@other_params$refuge_params$method, "competitive")
    expect_false(isTRUE(params@other_params$refuge_params$degrade))

    result <- reefDegrade(
        params,
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other,
        t = 1
    )
    expect_identical(result, params@other_params$refuge_params$method_params$refuge_density)
})

test_that("reefDegrade returns NULL for a non-competitive method", {
    data(caribbean_3_model)
    data(tuning_profile)
    params <- newRefuge(caribbean_3_model,
        new_method = "binned",
        new_method_params = tuning_profile
    )
    result <- reefDegrade(
        params,
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other,
        t = 1
    )
    expect_null(result)
})
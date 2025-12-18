test_that("reefDegrade returns NULL when degrade is FALSE", {
    data(caribbean_3_model)
    params <- caribbean_3_model

    result <- reefDegrade(
        params,
        n = params@initial_n,
        n_pp = params@initial_n_pp,
        n_other = params@initial_n_other,
        t = 1
    )
    expect_null(result)
})
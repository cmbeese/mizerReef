test_that("getSenMort forwards to reefSenMort with the given arguments", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n * 0.5
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other

    expect_equal(
        getSenMort(params, n = n, n_pp = n_pp, n_other = n_other, t = 3),
        reefSenMort(params, n = n, n_pp = n_pp, n_other = n_other, t = 3)
    )
})

test_that("getSenMort defaults to the params' own initial state", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_equal(
        getSenMort(params),
        reefSenMort(params,
            n = params@initial_n, n_pp = params@initial_n_pp,
            n_other = params@initial_n_other, t = 0
        )
    )
})

test_that("getSenMort rejects an n array of the wrong dimensions", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    bad_n <- params@initial_n[1, , drop = FALSE]
    expect_error(getSenMort(params, n = bad_n))
})

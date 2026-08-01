test_that("detritus_lifetime is exactly the reciprocal of detritus_consumption", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_equal(detritus_lifetime(params), 1 / detritus_consumption(params))
})

test_that("assigning a new detritus_lifetime rescales detritus so the new lifetime is achieved", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    target <- 2 * detritus_lifetime(params)

    detritus_lifetime(params) <- target
    expect_equal(detritus_lifetime(params), target)
})

test_that("assigning a new detritus_lifetime is equivalent to calling rescale_detritus() with the implied factor", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    target <- 3 * detritus_lifetime(params)
    implied_factor <- target / detritus_lifetime(params)
    expected <- rescale_detritus(params, implied_factor)

    detritus_lifetime(params) <- target
    expect_equal(params, expected)
})

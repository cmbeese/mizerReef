test_that("scaleReefModel runs without error", {
    params <- newReefParams()
    expect_error(scaleReefModel(params, 1), NA)
})

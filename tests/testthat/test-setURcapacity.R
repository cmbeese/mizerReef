test_that("setURcapacity runs without error", {
    params <- newReefParams()
    expect_error(setURcapacity(params), NA)
})

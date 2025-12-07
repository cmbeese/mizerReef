test_that("removeSpecies runs without error", {
    params <- newReefParams()
    expect_error(removeSpecies(params, 1), NA)
})

test_that("removeSpecies runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(removeSpecies(params, 1), NA)
})

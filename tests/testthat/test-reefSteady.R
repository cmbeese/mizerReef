test_that("reefSteady runs without error", {
    data(caribbean_3_model)
    expect_error(reefSteady(caribbean_3_model), NA)
})

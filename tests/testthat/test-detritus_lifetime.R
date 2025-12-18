test_that("detritus_lifetime runs without error", {
    data(caribbean_3_model)
    expect_error(detritus_lifetime(caribbean_3_model), NA)
})

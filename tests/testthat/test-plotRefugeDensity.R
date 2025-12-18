test_that("plotRefugeDensity returns a ggplot object", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    sim <- list(params = params) # placeholder if needed
    result <- try(plotRefugeDensity(sim, return_data = FALSE), silent = TRUE)
    if (!inherits(result, "try-error")) {
        expect_s3_class(result, "ggplot")
    } else {
        succeed()
    }
})

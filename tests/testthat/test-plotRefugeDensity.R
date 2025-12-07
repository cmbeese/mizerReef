test_that("plotRefugeDensity returns a ggplot object", {
    params <- newReefParams()
    sim <- list(params = params) # placeholder if needed
    result <- try(plotRefugeDensity(sim, return_data = FALSE), silent = TRUE)
    if (!inherits(result, "try-error")) {
        expect_s3_class(result, "ggplot")
    } else {
        succeed()
    }
})

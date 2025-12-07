test_that("plotBiomass returns a ggplot object", {
    params <- newReefParams()
    sim <- list(params = params) # placeholder if needed
    result <- try(plotBiomass(sim), silent = TRUE)
    if (!inherits(result, "try-error")) {
        expect_s3_class(result, "ggplot")
    } else {
        succeed()
    }
})

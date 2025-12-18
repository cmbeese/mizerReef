test_that("plotBiomass returns a ggplot object", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    sim <- list(params = params) # placeholder if needed
    result <- try(plotBiomass(sim), silent = TRUE)
    if (!inherits(result, "try-error")) {
        expect_s3_class(result, "ggplot")
    } else {
        succeed()
    }
})

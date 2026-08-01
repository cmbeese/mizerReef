test_that("scaleReefBackground matches scaleReefAbundance() followed by mizer::scaleModel()", {
    # scaleReefBackground() deliberately uses mizer's own scaleModel(), not
    # this package's scaleReefModel() -- it only rescales the background/
    # abundance, not algae or detritus, so confirm that's really what's
    # wired up.
    data(caribbean_3_model)
    params <- caribbean_3_model
    expected <- mizer::scaleModel(scaleReefAbundance(params, factor = 2), factor = 1 / 2)

    result <- scaleReefBackground(params, factor = 2)
    expect_equal(result, expected)
})

test_that("scaleReefBackground leaves algae/detritus parameters untouched", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- scaleReefBackground(params, factor = 2)

    expect_equal(result@other_params$algae$rho, params@other_params$algae$rho)
    expect_equal(result@other_params$detritus$rho, params@other_params$detritus$rho)
})

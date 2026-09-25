test_that("scaleReefBackground matches scaleReefAbundance() followed by scaleReefModel()", {
    # scaleReefBackground() calls mizer's scaleModel(), which on a mizerReef
    # model dispatches to scaleModel.mizerReef() and so scales algae and
    # detritus exactly as scaleReefModel() does.
    data(caribbean_3_model)
    params <- caribbean_3_model
    expected <- scaleReefModel(scaleReefAbundance(params, factor = 2), factor = 1 / 2)

    result <- scaleReefBackground(params, factor = 2)
    expect_equal(result, expected)
})

test_that("scaleReefBackground leaves the algae and detritus encounter rates unchanged", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- scaleReefBackground(params, factor = 2)

    for (component in c("algae", "detritus")) {
        expect_equal(
            encounter_contribution(result, result@initial_n_other, component),
            encounter_contribution(params, params@initial_n_other, component)
        )
    }
})

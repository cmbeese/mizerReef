test_that("getProductivity(params) matches an independently computed per-species sum", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    energy <- getEGrowthTime(params)
    size_range <- mizer::get_size_range_array(params,
        min_l = 7, max_l = max(params@species_params$l_max)
    )
    n <- params@initial_n

    expected <- vapply(seq_len(nrow(n)), function(sp) {
        sum(energy[sp, ] * n[sp, ] * size_range[sp, ] * params@dw)
    }, numeric(1))
    names(expected) <- params@species_params$species

    expect_equal(as.numeric(getProductivity(params)), as.numeric(expected))
})

test_that("getProductivity defaults min_fishing_l to 7cm and max_fishing_l to the largest species' l_max", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    explicit <- getProductivity(params,
        min_fishing_l = 7, max_fishing_l = max(params@species_params$l_max)
    )
    expect_equal(getProductivity(params), explicit)
})

test_that("getProductivity with a narrower size range gives less than or equal productivity", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    full <- getProductivity(params)
    narrow <- getProductivity(params, min_fishing_l = 7, max_fishing_l = 10)

    expect_true(all(as.numeric(narrow) <= as.numeric(full) + 1e-8))
    expect_true(any(as.numeric(narrow) < as.numeric(full)))
})

test_that("getProductivity with include_repro = TRUE uses getEReproAndGrowth instead of getEGrowthTime", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    energy <- mizer::getEReproAndGrowth(params)
    size_range <- mizer::get_size_range_array(params,
        min_l = 7, max_l = max(params@species_params$l_max)
    )
    n <- params@initial_n
    expected <- vapply(seq_len(nrow(n)), function(sp) {
        sum(energy[sp, ] * n[sp, ] * size_range[sp, ] * params@dw)
    }, numeric(1))

    result <- getProductivity(params, include_repro = TRUE)
    expect_equal(as.numeric(result), expected)
    expect_false(isTRUE(all.equal(as.numeric(result), as.numeric(getProductivity(params)))))
})

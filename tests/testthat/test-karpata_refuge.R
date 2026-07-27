test_that("karpata_refuge is a valid competitive-method method_params data frame", {
    data(karpata_refuge)
    expect_s3_class(karpata_refuge, "data.frame")
    expect_true(all(c("start_L", "end_L", "refuge_density") %in% names(karpata_refuge)))
    expect_true(all(karpata_refuge$start_L < karpata_refuge$end_L))
    expect_true(all(karpata_refuge$refuge_density >= 0))
    expect_false(anyNA(karpata_refuge))
})

test_that("karpata_refuge works as setRefuge()/getRefuge() input for the competitive method", {
    data(caribbean_10_species)
    data(caribbean_10_interaction)
    data(karpata_refuge)
    base <- suppressMessages(newMultispeciesParams(
        species_params = caribbean_10_species, interaction = caribbean_10_interaction,
        min_w_pp = NA, w_pp_cutoff = 1, n = 0.75, p = 0.75
    ))
    base@other_params$refuge_params <- list()
    result <- getRefuge(setRefuge(base, method = "competitive", method_params = karpata_refuge))
    expect_equal(
        result@other_params$refuge_params$method_params,
        karpata_refuge,
        ignore_attr = TRUE
    )
})

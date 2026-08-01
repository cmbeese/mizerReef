test_that("getRefuge computes the documented sigmoidal formula with dummy fish bins", {
    # ref = prop_protect / (1 + exp(slope*(w - W_refuge))), capped at
    # max_protect, zeroed below w_settle, and zeroed for non-refuge-users --
    # hand-computed here directly from w/a_bar/b_bar, not by re-deriving
    # getRefuge()'s own code.
    data(caribbean_3_model)
    params <- caribbean_3_model
    method_params <- data.frame(L_refuge = 5, prop_protect = 0.6)
    result <- getRefuge(setRefuge(params, method = "sigmoidal", method_params = method_params))

    w <- params@w
    a_bar <- 0.025
    b_bar <- 3
    w_settle <- 0.1
    max_protect <- 0.98
    slope <- 100
    W_refuge <- a_bar * method_params$L_refuge^b_bar
    denom <- 1 + exp(slope * (w - W_refuge))
    ref <- ifelse(w > w_settle, method_params$prop_protect / denom, 0)
    ref[ref > max_protect] <- max_protect
    refuge_user <- params@species_params$refuge_user
    expected <- matrix(0, nrow = 3, ncol = length(w))
    for (i in 1:3) expected[i, ] <- refuge_user[i] * ref

    expect_equal(unname(result@other_params$refuge_params$refuge), unname(expected))
})

test_that("getRefuge computes the documented sigmoidal formula with species-specific bins", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    method_params <- data.frame(L_refuge = 5, prop_protect = 0.6)
    result <- getRefuge(
        setRefuge(params,
            method = "sigmoidal", method_params = method_params,
            use_dummy_fish_bins = FALSE
        ),
        use_dummy_fish_bins = FALSE
    )

    w <- params@w
    sp <- params@species_params
    w_settle <- 0.1
    max_protect <- 0.98
    slope <- 100
    expected <- matrix(0, nrow = 3, ncol = length(w))
    for (i in 1:3) {
        W_refuge_i <- sp$a[i] * method_params$L_refuge^sp$b[i]
        denom <- 1 + exp(slope * (w - W_refuge_i))
        ref <- ifelse(w > w_settle, method_params$prop_protect / denom, 0)
        ref[ref > max_protect] <- max_protect
        expected[i, ] <- sp$refuge_user[i] * ref
    }

    expect_equal(unname(result@other_params$refuge_params$refuge), unname(expected))
})

test_that("getRefuge computes the documented binned formula with dummy fish bins", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    method_params <- data.frame(
        start_L = c(0, 10, 20), end_L = c(10, 20, 30), prop_protect = c(0.8, 0.5, 0.2)
    )
    result <- getRefuge(setRefuge(params, method = "binned", method_params = method_params))

    w <- params@w
    a_bar <- 0.025
    b_bar <- 3
    w_settle <- 0.1
    max_protect <- 0.98
    ref <- rep(0, length(w))
    for (k in 1:3) {
        start_w <- max(a_bar * method_params$start_L[k]^b_bar, w_settle)
        end_w <- a_bar * method_params$end_L[k]^b_bar
        ref[which(w >= start_w & w <= end_w)] <- method_params$prop_protect[k]
    }
    refuge_user <- params@species_params$refuge_user
    expected <- matrix(0, nrow = 3, ncol = length(w))
    for (i in 1:3) expected[i, ] <- refuge_user[i] * pmin(ref, max_protect)

    expect_equal(unname(result@other_params$refuge_params$refuge), unname(expected))
})

test_that("getRefuge computes the documented binned formula with species-specific bins", {
    # caribbean_3_model's 3 species all share a = 0.025, b = 3, so it can't
    # distinguish species-specific bins from dummy bins. caribbean_10_species
    # has genuinely different a/b per species, which is needed to catch a
    # real bug this test originally found: the species-specific branch used
    # to share a single refuge-proportion vector across every species/bin
    # combination (so every species ended up with an identical profile,
    # defeating the whole point of use_dummy_fish_bins = FALSE) and also
    # crashed building refuge_lengths whenever there was more than one
    # species. Both are fixed; this checks the corrected numeric output.
    data(caribbean_10_species)
    data(caribbean_10_interaction)
    sp <- caribbean_10_species
    base <- suppressMessages(newMultispeciesParams(
        species_params = sp, interaction = caribbean_10_interaction,
        min_w_pp = NA, w_pp_cutoff = 1, n = 0.75, p = 0.75
    ))
    base@other_params$refuge_params <- list()
    method_params <- data.frame(
        start_L = c(0, 10, 20), end_L = c(10, 20, 30), prop_protect = c(0.8, 0.5, 0.2)
    )
    result <- getRefuge(
        setRefuge(base, method = "binned", method_params = method_params, use_dummy_fish_bins = FALSE),
        use_dummy_fish_bins = FALSE
    )
    refuge_actual <- result@other_params$refuge_params$refuge

    w <- base@w
    w_settle <- 0.1
    max_protect <- 0.98
    no_sp <- nrow(sp)
    expected <- matrix(0, nrow = no_sp, ncol = length(w))
    for (i in seq_len(no_sp)) {
        ref_i <- rep(0, length(w))
        for (k in 1:3) {
            start_w <- max(sp$a[i] * method_params$start_L[k]^sp$b[i], w_settle)
            end_w <- sp$a[i] * method_params$end_L[k]^sp$b[i]
            ref_i[which(w >= start_w & w <= end_w)] <- method_params$prop_protect[k]
        }
        expected[i, ] <- sp$refuge_user[i] * pmin(ref_i, max_protect)
    }

    expect_equal(unname(refuge_actual), unname(expected))
    # The bug this guards against produced identical rows for every species;
    # confirm the fixed output is genuinely species-specific.
    expect_false(isTRUE(all.equal(refuge_actual[1, ], refuge_actual[3, ])))
})

test_that("getRefuge computes competitive method bin.id with species-specific bins", {
    data(caribbean_10_species)
    data(caribbean_10_interaction)
    sp <- caribbean_10_species
    base <- suppressMessages(newMultispeciesParams(
        species_params = sp, interaction = caribbean_10_interaction,
        min_w_pp = NA, w_pp_cutoff = 1, n = 0.75, p = 0.75
    ))
    base@other_params$refuge_params <- list()
    method_params <- data.frame(
        start_L = c(0, 10, 20), end_L = c(10, 20, 30), refuge_density = c(1, 0.5, 0.1)
    )
    result <- getRefuge(
        setRefuge(base, method = "competitive", method_params = method_params, use_dummy_fish_bins = FALSE),
        use_dummy_fish_bins = FALSE
    )
    bin_id <- result@other_params$refuge_params$bin.id

    w <- base@w
    w_settle <- 0.1
    for (i in seq_len(nrow(sp))) {
        for (k in 1:3) {
            key <- paste0("sp", i, "_bin", k)
            start_w <- max(sp$a[i] * method_params$start_L[k]^sp$b[i], w_settle)
            end_w <- sp$a[i] * method_params$end_L[k]^sp$b[i]
            expected_idx <- which(w >= start_w & w <= end_w)
            expect_equal(bin_id[[key]], expected_idx, info = key)
        }
    }
})

test_that("getRefuge computes competitive method bin.id as the weight range for each length bin", {
    # caribbean_3_model already uses the competitive method with a 10-bin
    # method_params -- rebuild bin.id independently from start_L/end_L/a_bar/b_bar.
    data(caribbean_3_model)
    params <- caribbean_3_model
    method_params <- params@other_params$refuge_params$method_params
    result <- getRefuge(params)
    bin_id <- result@other_params$refuge_params$bin.id

    w <- params@w
    a_bar <- 0.025
    b_bar <- 3
    w_settle <- 0.1
    for (k in seq_len(nrow(method_params))) {
        start_w <- max(a_bar * method_params$start_L[k]^b_bar, w_settle)
        end_w <- a_bar * method_params$end_L[k]^b_bar
        expected_idx <- which(w >= start_w & w <= end_w)
        expect_equal(bin_id[[k]], expected_idx, info = paste("bin", k))
    }
})

test_that("getRefuge produces an all-zero refuge matrix for the noncomplex method", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- getRefuge(setRefuge(params, method = "noncomplex"))
    refuge <- result@other_params$refuge_params$refuge
    expect_equal(dim(refuge), c(3L, length(params@w)))
    expect_true(all(refuge == 0))
})

test_that("getRefuge uses use_dummy_fish_bins stored in params when not passed explicitly", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    method_params <- data.frame(L_refuge = 5, prop_protect = 0.6)

    stored_false <- setRefuge(params,
        method = "sigmoidal", method_params = method_params, use_dummy_fish_bins = FALSE
    )
    result_default <- getRefuge(stored_false, use_dummy_fish_bins = NULL)
    result_explicit <- getRefuge(stored_false, use_dummy_fish_bins = FALSE)

    expect_equal(
        result_default@other_params$refuge_params$refuge,
        result_explicit@other_params$refuge_params$refuge
    )
})

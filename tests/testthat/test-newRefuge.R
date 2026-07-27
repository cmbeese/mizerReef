test_that("newRefuge with sigmoidal method matches setRefuge + getRefuge called directly", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    method_params <- data.frame(L_refuge = 0.5, prop_protect = 0.3)

    via_new_refuge <- newRefuge(params,
        new_method = "sigmoidal", new_method_params = method_params
    )
    via_direct_calls <- getRefuge(setRefuge(params,
        method = "sigmoidal", method_params = method_params
    ))

    expect_equal(
        via_new_refuge@other_params$refuge_params$refuge,
        via_direct_calls@other_params$refuge_params$refuge
    )
    expect_equal(
        via_new_refuge@other_params$refuge_params$refuge_lengths,
        via_direct_calls@other_params$refuge_params$refuge_lengths
    )
})

test_that("newRefuge with sigmoidal method and list params matches the data frame version", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    mp_df <- data.frame(L_refuge = 0.5, prop_protect = 0.3)
    mp_list <- list(L_refuge = 0.5, prop_protect = 0.3)

    from_df <- newRefuge(params, new_method = "sigmoidal", new_method_params = mp_df)
    from_list <- newRefuge(params, new_method = "sigmoidal", new_method_params = mp_list)

    expect_equal(
        from_df@other_params$refuge_params$refuge,
        from_list@other_params$refuge_params$refuge
    )
})

test_that("newRefuge sets the new_refuge flag as requested", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    method_params <- data.frame(L_refuge = 0.5, prop_protect = 0.3)

    result_true <- newRefuge(params,
        new_refuge = TRUE,
        new_method = "sigmoidal", new_method_params = method_params
    )
    result_false <- newRefuge(params,
        new_refuge = FALSE,
        new_method = "sigmoidal", new_method_params = method_params
    )
    expect_true(result_true@other_params$new_refuge)
    expect_false(result_false@other_params$new_refuge)
})

test_that("newRefuge scale_bin scales refuge_density for the competitive method", {
    # caribbean_3_model uses the competitive method with a 10-row method_params
    data(caribbean_3_model)
    params <- caribbean_3_model
    before_rd <- params@other_params$refuge_params$method_params$refuge_density

    scaled_scalar <- suppressWarnings(newRefuge(params, scale_bin = 2))
    expect_equal(
        scaled_scalar@other_params$refuge_params$method_params$refuge_density,
        before_rd * 2
    )

    scale_vec <- seq(1, 2, length.out = length(before_rd))
    scaled_vector <- suppressWarnings(newRefuge(params, scale_bin = scale_vec))
    expect_equal(
        scaled_vector@other_params$refuge_params$method_params$refuge_density,
        before_rd * scale_vec
    )
})

test_that("newRefuge scale_bin scales prop_protect for the binned method", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    method_params <- data.frame(
        start_L = 1:10, end_L = 2:11, prop_protect = rep(0.5, 10)
    )
    binned <- newRefuge(params, new_method = "binned", new_method_params = method_params)
    before_pp <- binned@other_params$refuge_params$method_params$prop_protect

    scaled_scalar <- suppressWarnings(newRefuge(binned, scale_bin = 2))
    expect_equal(
        scaled_scalar@other_params$refuge_params$method_params$prop_protect,
        before_pp * 2
    )

    scale_vec <- seq(1, 1.5, length.out = length(before_pp))
    scaled_vector <- suppressWarnings(newRefuge(binned, scale_bin = scale_vec))
    expect_equal(
        scaled_vector@other_params$refuge_params$method_params$prop_protect,
        before_pp * scale_vec
    )
})

test_that("newRefuge errors when no input is provided", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(newRefuge(params), "At least one input must be provided")
})

test_that("newRefuge errors when switching method without new_method_params", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(
        newRefuge(params, new_method = "sigmoidal"),
        "must provide a new method_params data frame to"
    )
})

test_that("newRefuge errors when scale_bin length matches neither 1 nor the bin count", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    no_bins <- nrow(params@other_params$refuge_params$method_params)
    expect_gt(no_bins, 1)
    expect_error(
        suppressWarnings(newRefuge(params, scale_bin = c(1, 2, 3))),
        "scale_bin must have length 1 or have"
    )
})

test_that("newRefuge errors when scale_bin is negative", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    no_bins <- nrow(params@other_params$refuge_params$method_params)
    expect_error(
        suppressWarnings(newRefuge(params, scale_bin = c(rep(1, no_bins - 1), -1))),
        "scale_bin must be non-negative"
    )
})

test_that("newRefuge warns and reuses the stored method when new_method is not given", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_warning(
        result <- newRefuge(params, scale_bin = 2),
        "did not provide a new method"
    )
    expect_identical(
        result@other_params$refuge_params$method,
        params@other_params$refuge_params$method
    )
})

test_that("setDegradation runs with a real trajectory and enables degradation", {
    data(caribbean_3_model)
    data(rubble_scale)
    params <- caribbean_3_model
    params <- setDegradation(params,
        trajectory = "rubble", deg_scale = rubble_scale,
        bleach_time = 2, degrade = TRUE
    )
    expect_true(params@other_params$refuge_params$degrade)
    expect_identical(params@other_params$refuge_params$trajectory, "rubble")

    # Exercise the actual degradation time-series logic in reefDegrade():
    # refuge density before the bleach event should differ from after it.
    rd_before <- reefDegrade(params,
        n = params@initial_n, n_pp = params@initial_n_pp,
        n_other = params@initial_n_other, t = 1
    )
    rd_after <- reefDegrade(params,
        n = params@initial_n, n_pp = params@initial_n_pp,
        n_other = params@initial_n_other, t = 3
    )
    expect_false(isTRUE(all.equal(rd_before, rd_after)))
})

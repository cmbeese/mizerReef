test_that("rubble_scale/algae_scale/recovery_scale are valid, consistently-shaped degradation trajectory matrices", {
    data(rubble_scale)
    data(algae_scale)
    data(recovery_scale)

    for (traj in list(rubble = rubble_scale, algae = algae_scale, recovery = recovery_scale)) {
        expect_true(is.matrix(traj))
        expect_true(is.numeric(traj))
        expect_false(anyNA(traj))
        expect_true(all(traj >= 0))
        expect_equal(dim(traj), c(10L, 15L))
    }

    expect_identical(rownames(rubble_scale), rownames(algae_scale))
    expect_identical(rownames(rubble_scale), rownames(recovery_scale))
    expect_identical(colnames(rubble_scale), colnames(algae_scale))
    expect_identical(colnames(rubble_scale), colnames(recovery_scale))
})

test_that("degradation trajectories work as setDegradation()/reefDegrade() input", {
    data(caribbean_3_model)
    data(rubble_scale)
    data(algae_scale)
    data(recovery_scale)
    params0 <- caribbean_3_model
    n <- params0@initial_n
    n_pp <- params0@initial_n_pp
    n_other <- params0@initial_n_other

    for (traj in list(rubble_scale, algae_scale, recovery_scale)) {
        params <- setDegradation(params0, deg_scale = traj, bleach_time = 2, degrade = TRUE)
        result <- reefDegrade(params, n = n, n_pp = n_pp, n_other = n_other, t = 2)
        expect_false(anyNA(result))
        expect_true(all(result >= 0))
    }
})

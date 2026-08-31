# mizer 3.3 routes what a model-building function tells the user through
# signal_info()/with_info_level(), controlled by an `info_level` argument.
# These tests check that mizerReef's own reports go through that mechanism
# rather than being bare message()/warning() calls the user cannot turn down.

warnings_from <- function(expr) {
    w <- character()
    withCallingHandlers(
        force(expr),
        warning = function(cnd) {
            w <<- c(w, conditionMessage(cnd))
            invokeRestart("muffleWarning")
        },
        message = function(cnd) invokeRestart("muffleMessage")
    )
    w
}

# caribbean_3_model is now well-balanced at steady state (its detritus
# scale/lifetime having just been tuned to a literature target), so it no
# longer incidentally triggers tuneUR()/tuneUR_cc()'s negative-external-flux
# case on its own the way it used to. Force a detritus consumption deficit
# deliberately -- `rho` is cached on `other_params$detritus$rho` (an outer
# product with w^n, not recomputed from `species_params$rho_detritus` on
# the fly), so both need scaling down together for `getDetritusConsumption()`
# to see it.
force_detritus_deficit <- function(params) {
    params@other_params$detritus$rho <- params@other_params$detritus$rho * 0.01
    params@species_params$rho_detritus <- params@species_params$rho_detritus * 0.01
    params
}

test_that("setAlgaeParams reports the interaction_algae default through info_level", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@species_params$interaction_algae <- NULL

    expect_match(warnings_from(setAlgaeParams(params)), "interaction_algae",
                 all = FALSE)
    expect_length(warnings_from(setAlgaeParams(params, info_level = 0)), 0)
})

test_that("setDetritusParams reports the interaction_detritus default through info_level", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@species_params$interaction_detritus <- NULL

    expect_match(warnings_from(setDetritusParams(params)), "interaction_detritus",
                 all = FALSE)
    expect_length(warnings_from(setDetritusParams(params, info_level = 0)), 0)
})

test_that("tuneUR and tuneUR_cc report the negative external detritus flux through info_level", {
    data(caribbean_3_model)
    params <- force_detritus_deficit(caribbean_3_model)

    expect_match(warnings_from(tuneUR(params)),
                 "flux of external detritus is negative", all = FALSE)
    expect_length(warnings_from(tuneUR(params, info_level = 0)), 0)

    cc_params <- suppressWarnings(setURcapacity(params, cap = 1.5))
    expect_length(warnings_from(tuneUR_cc(cc_params, info_level = 0)), 0)
})

test_that("reefSteady carries info_level down to the algae/detritus tuning", {
    # The reports raised by tuneUR()/tuneUR_cc() are collected by the handler
    # reefSteady() installs, so a caller can quiet the whole call.
    data(caribbean_3_model)
    params <- force_detritus_deficit(caribbean_3_model)
    expect_match(warnings_from(reefSteady(params, progress_bar = FALSE)),
                 "flux of external detritus is negative", all = FALSE)
    expect_length(
        warnings_from(reefSteady(params, progress_bar = FALSE,
                                 info_level = 0)),
        0
    )
})

test_that("newReefParams collects the reports of the setters it calls", {
    # newReefParams() installs one handler for the whole construction and
    # also forwards info_level explicitly to setRefuge()/setAlgaeParams()/
    # setDetritusParams() -- see the regression test below for why forwarding
    # matters even though nested with_info_level() calls "just work" in the
    # common case.
    data(caribbean_3_species)
    data(caribbean_3_interaction)
    data(tuning_profile)
    build <- function(...) {
        newReefParams(species_params = caribbean_3_species,
                      interaction = caribbean_3_interaction,
                      method = "binned", method_params = tuning_profile, ...)
    }
    expect_length(warnings_from(build(info_level = 0)), 0)
    expect_s4_class(suppressMessages(suppressWarnings(build(info_level = 0))),
                    "mizerReef")
})

test_that("an explicit info_level overrides a differing global mizer_info_level option", {
    # mizer::with_info_level()'s nesting only "just works" without forwarding
    # when an inner call's own resolved info_level happens to agree with the
    # outer one -- e.g. both left at the default, which reads the same
    # global option either way. If a global `options(mizer_info_level = 0)`
    # differs from what the outer call was explicitly given, an unforwarded
    # inner call resolves its *own* default independently, gets 0, and takes
    # with_info_level()'s documented "silence is the exception" path:
    # unconditionally muffling its reports regardless of what the outer,
    # explicit info_level asked for. reefSteady()/newReefParams()/newRefuge()
    # forward info_level to their inner setters/tuners specifically to avoid
    # this.
    data(caribbean_3_model)
    withr::local_options(mizer_info_level = 0)

    expect_match(
        warnings_from(reefSteady(force_detritus_deficit(caribbean_3_model),
                                 progress_bar = FALSE, info_level = 3)),
        "flux of external detritus is negative", all = FALSE
    )

    params <- caribbean_3_model
    params@species_params$interaction_algae <- NULL
    expect_match(
        warnings_from(setAlgaeParams(params, info_level = 3)),
        "interaction_algae", all = FALSE
    )
})

test_that("the reports keep the wording they had as bare warnings", {
    # Converting to signal_info() must not change what the user reads, so that
    # existing scripts matching on these messages keep working.
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@species_params$refuge_user <- NULL

    w <- warnings_from(setRefuge(params,
                                 method = params@other_params$refuge_params$method,
                                 method_params = params@other_params$refuge_params$method_params))
    expect_match(w, "You have not provided values for refuge_user", all = FALSE)
})

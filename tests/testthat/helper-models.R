# caribbean_3_model with linear algae and detritus dynamics (no carrying
# capacity), with algae and detritus balanced for those dynamics. For tests of
# the linear dynamics and tuneUR(): the bundled model gives algae a carrying
# capacity.
linear_caribbean_3 <- function() {
    data("caribbean_3_model", package = "mizerReef", envir = environment())
    params <- caribbean_3_model
    params@other_params$use_UR_cc <- FALSE
    params@other_dynamics$algae <- "algae_dynamics"
    params@other_dynamics$detritus <- "detritus_dynamics"
    tuneUR(params, info_level = 0)
}

# The bundled models were calibrated with the senescence form mizerReef used
# before 2.1.0 (sen_form = "log_ratio"); saved models now use the power form.
# Until the models are recalibrated, tests of their steady state put back
# the form they were calibrated with.
with_calibration_senescence <- function(params) {
    withr::with_options(
        list(lifecycle_verbosity = "quiet"),
        setExtMortParams(params, params@other_params$ext_mort_params,
                         sen_form = "log_ratio"))
}

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

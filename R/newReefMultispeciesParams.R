#' Set up parameters for a mizerReef model
#'
#' @return An object of type \linkS4class{MizerParams}
#'
#' @export
#' @family functions for setting up models
newReefMultispeciesParams <- function(species_params,
                                      interaction = NULL,
                                      no_w = 100,
                                      min_w = 0.001,
                                      max_w = NA,
                                      min_w_pp = NA,
                                      # setPredKernel()
                                      pred_kernel = NULL,
                                      # setSearchVolume()
                                      search_vol = NULL,
                                      # setMaxIntakeRate()
                                      intake_max = NULL,
                                      # setMetabolicRate()
                                      metab = NULL,
                                      p = 0.7,
                                      # setExtMort
                                      ext_mort = NULL,
                                      z0pre = 0.6,
                                      z0exp = n - 1,
                                      # setReproduction
                                      maturity = NULL,
                                      repro_prop = NULL,
                                      RDD = "BevertonHoltRDD",
                                      # setResource
                                      kappa = 1e11,
                                      n = 2 / 3,
                                      resource_rate = 10,
                                      resource_capacity = kappa,
                                      lambda = 2.05,
                                      w_pp_cutoff = 10,
                                      resource_dynamics = "resource_semichemostat",
                                      gear_params = NULL,
                                      selectivity = NULL,
                                      catchability = NULL,
                                      initial_effort = NULL,
                                      info_level = 3,
                                      z0 = lifecycle::deprecated(),
                                      r_pp = lifecycle::deprecated()) {

    if (lifecycle::is_present(r_pp)) {
        lifecycle::deprecate_warn("1.0.0", "newMultispeciesParams(r_pp)",
                                  "newMultispeciesParams(resource_rate)")
        resource_rate <- r_pp
    }
    if (lifecycle::is_present(z0)) {
        lifecycle::deprecate_warn("2.2.3", "newMultispeciesParams(z0)",
                                  "newMultispeciesParams(ext_mort)")
        ext_mort <- z0
    }

    # Define a signal handler that collects the information signals
    # into the `infos` list.
    infos <- list()
    collect_info <- function(cnd) {
        if (cnd$level <= info_level) {
            infos[[cnd$var]] <<- cnd$message
        }
    }
    # Register this signal handler
    withCallingHandlers(
        info_about_default = collect_info, {
            no_sp <- nrow(species_params)
            species_params <- validSpeciesParams(species_params)
            gear_params <- validGearParams(gear_params, species_params)

            ## Create MizerReefParams object ----
            params <- emptyReefParams(species_params,
                                      gear_params,
                                      no_w = no_w,
                                      min_w = min_w,
                                      max_w = max_w,
                                      min_w_pp = min_w_pp)

            # Fill the slots ----
            params <- params %>%
                set_species_param_default("n", n) %>%
                set_species_param_default("p", p)
            params <- set_species_param_default(
                params, "q", lambda - 2 + params@species_params[["n"]])
            if (is.null(interaction)) {
                interaction <- matrix(1, nrow = no_sp, ncol = no_sp)
            }

            params@initial_n_pp[] <- kappa * params@w_full ^ (-lambda)
            params@initial_n_pp[params@w_full >= w_pp_cutoff] <- 0
            params@resource_params$kappa <- kappa
            params@resource_params$lambda <- lambda
            params@resource_params$w_pp_cutoff <- w_pp_cutoff

            params <- params  %>%
                setParams(
                    # setInteraction
                    interaction = interaction,
                    # setPredKernel()
                    pred_kernel = pred_kernel,
                    # setSearchVolume()
                    search_vol = search_vol,
                    # setMaxIntakeRate()
                    intake_max = intake_max,
                    # setMetabolicRate()
                    metab = metab,
                    # setExtMort
                    ext_mort = ext_mort,
                    z0pre = z0pre,
                    z0exp = z0exp,
                    # setReproduction
                    maturity = maturity,
                    repro_prop = repro_prop,
                    RDD = RDD,
                    # setFishing
                    gear_params = gear_params,
                    selectivity = selectivity,
                    catchability = catchability,
                    initial_effort = initial_effort) %>%
                setResource(
                    # setResource
                    resource_rate = resource_rate,
                    resource_capacity = resource_capacity,
                    resource_dynamics = resource_dynamics,
                    lambda = lambda,
                    n = n,
                    w_pp_cutoff = w_pp_cutoff,
                    balance = FALSE)

            params@initial_n <- mizer::get_initial_n(params)
            params@A <- rep(1, nrow(species_params))
        })
    if (length(infos) > 0) {
        message(paste(infos, collapse = "\n"))
    }
    return(params)
}


# helper function to calculate w_min_idx slot
get_w_min_idx <- function(species_params, w) {
    # Round down w_min to lie on grid points and store the indices of these
    # grid points in w_min_idx
    w_min_idx <- as.vector(suppressWarnings(
        tapply(species_params$w_min, seq_len(nrow(species_params)),
               function(w_min, wx) max(which(wx <= w_min)), wx = w)))
    # Due to rounding errors this might happen:
    w_min_idx[w_min_idx == -Inf] <- 1
    names(w_min_idx) <- as.character(species_params$species)
    w_min_idx
}

# environment(newReefMultispeciesParams) <- asNamespace("mizer")
# utils::assignInNamespace("newMultispeciesParams", newReefMultispeciesParams, ns = "mizer")

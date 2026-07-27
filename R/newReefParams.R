#' Set up parameters for a mizerReef model
#'
#' Sets up a multi-species size spectrum model with additional unstructured
#' resource components, senescence mortality, and predation refuge.
#'
#' @inheritSection setAlgaeParams Algae as an unstructured resource
#' @inheritSection setDetritusParams Detritus as an unstructured resource
#' @inheritSection setExtMortParams Senescence mortality
#' @inheritSection setRefuge Setting the refuge profile
#'
#' @export
#' @inheritParams setAlgaeParams
#' @inheritParams setDetritusParams
#' @inheritParams setExtMortParams
#' @inheritParams setRefuge
#' @inheritParams setDegradation
#'
#' @param species_params A species parameter data frame containing at
#'                       least the name of each species, their
#'                       observed abundances, and their maximum size.
#'
#' @param interaction The species specific interaction matrix, \eqn{\theta_{ij}}
#'
#' @param include_ext_mort A boolean value that indicates whether the user wants
#'                         to use default external mortality. Defaults to TRUE.
#'
#' @param include_sen_mort A boolean value that indicates whether the user wants
#'                         to use default senescence mortality. Defaults to TRUE.
#'
#' @param z0pre If `include_ext_mort`is FALSE, the external mortality rate for
#'              each species calculated as z0pre * w_max ^ z0exp. z0exp defaults
#'              to 1-n where n is the given allometric scaling exponent and
#'              z0pre defaults to 0.2.
#'
#' @param new_refuge Logical. If TRUE, indicates this refuge profile is being
#'                   set for use in a new simulation (prevents algae/detritus
#'                   from being re-tuned during [reefSteady()]). Default FALSE.
#'
#' @param degrade Logical. Whether to enable habitat degradation during
#'                projections. Default FALSE. See [setDegradation()] for
#'                details on degradation parameters.
#'
#' @param resource_color Character. Colour to use for the resource line in plots.
#'                       Default is "lightseagreen".
#'
#' @param initial_algae_growth Numeric. The initial growth rate of algae in
#'                             grams/m^2/year, passed on to [setAlgaeParams()]
#'                             as `algae_growth_initial`. Default is 2e3.
#'
#' @param min_w_pp Minimum size of plankton in grams
#'
#' @param w_pp_cutoff Maximum size of plankton in grams, default to 1 g
#'
#' @param n Allometric growth exponent (also used as metabolic exponent p)
#'
#' @param crit_feed Critical feeding level
#'
#' @param ... Extra parameters to be passed to [newMultispeciesParams()]
#'
#' @concept setup
#' @return An object of type [MizerParams]
newReefParams <- function( # Original mizer parameters
                          species_params, interaction = NULL,
                          crit_feed = NULL,
                          min_w_pp = NA, w_pp_cutoff = 1,
                          n = 0.75,
                          # Parameters for setting up refuge
                          new_refuge = FALSE,
                          method, method_params,
                          refuge_user = NULL, blocked_pred = NULL,
                          satiation = NULL,
                          a_bar = NULL, b_bar = NULL,
                          w_settle = NULL, max_protect = NULL,
                          tau = NULL,
                          use_dummy_fish_bins = TRUE,
                          # Parameters for setting up degradation
                          degrade = FALSE, bleach_time = 2,
                          trajectory = NULL, deg_scale = 1,
                          algae_boost = FALSE,
                          algae_growth_boost = c(1.11, 1.11, 1.11, 1.11),
                          algae_capacity_boost = c(2.0),
                          # Parameters for unstructured resources
                          UR_interaction = NULL,
                          use_UR_cc = FALSE,
                          initial_algae_growth = NULL,
                          algae_capacity = NULL,
                          detritus_capacity = NULL,
                          sen_decomp = NULL, ext_decomp = NULL,
                          external = NULL,
                          # Resource aesthetics
                          algae_colour = "darkseagreen3",
                          detritus_colour = "plum4",
                          resource_color = "lightseagreen",
                          # Parameters for external mortality
                          ext_mort_params = NULL,
                          include_ext_mort = TRUE,
                          include_sen_mort = TRUE,
                          z0pre = 0.2, ...) {
    ## Initialize model with newMultispeciesParams ----
    params <- newMultispeciesParams(
        species_params = species_params,
        interaction = interaction,
        min_w_pp = min_w_pp,
        w_pp_cutoff = w_pp_cutoff,
        n = n, p = n, ...
    )

    # Initialise other_params sub-lists for reef-specific state
    # (these replace the old custom S4 slots)
    if (is.null(params@other_params$refuge_params))  params@other_params$refuge_params  <- list()
    if (is.null(params@other_params$algae_params))   params@other_params$algae_params   <- list()
    if (is.null(params@other_params$detritus_params)) params@other_params$detritus_params <- list()

    # Change resource color
    params@linecolour["Resource"] <- resource_color
    params@linecolour["Fishing"] <- "blue"

    # Save include_sen_mort
    params@other_params$include_sen_mort <- include_sen_mort

    # Add parameters ----

    ### Refuge ----
    # Set new_refuge to FALSE for initial models
    # (this flag prevents recalculation of algae & detritus dynamics when changing refuge and running reefSteady)
    params@other_params$new_refuge <- FALSE

    # Add refuge parameters
    params <- setRefuge(
        params = params, method = method,
        method_params = method_params,
        refuge_user = refuge_user,
        blocked_pred = blocked_pred,
        satiation = satiation,
        a_bar = a_bar, b_bar = b_bar,
        w_settle = w_settle,
        max_protect = max_protect, tau = tau,
        use_dummy_fish_bins = use_dummy_fish_bins, ...
    )

    # Find initial refuge profiles
    params <- getRefuge(params, ...)

    ### Unstructured resources ----
    #### Algae ----
    params <- setAlgaeParams(
        params = params,
        algae_growth_initial = initial_algae_growth,
        algae_capacity = algae_capacity,
        UR_interaction = UR_interaction,
        use_UR_cc = use_UR_cc,
        algae_colour = algae_colour
    )

    #### Detritus ----
    params <- setDetritusParams(
        params = params,
        detritus_capacity = detritus_capacity,
        sen_decomp = sen_decomp,
        ext_decomp = ext_decomp,
        external = external,
        UR_interaction = UR_interaction,
        use_UR_cc = use_UR_cc,
        detritus_colour = detritus_colour
    )

    ### External mortality ----
    params <- setExtMortParams(
        params = params,
        ext_mort_params = ext_mort_params
    )

    ### Degradation ----
    if (isTRUE(degrade)) {
        params <- setDegradation(params,
            deg_scale = deg_scale,
            bleach_time = bleach_time,
            trajectory = trajectory,
            degrade = degrade,
            algae_boost = algae_boost,
            algae_growth_boost = algae_growth_boost,
            algae_capacity_boost = algae_capacity_boost
        )
    } else {
        params@other_params$refuge_params$degrade <- FALSE
    }

    # Algae & Detritus ----
    ### Calculate Rho ----
    # Determine the necessary detritus and algae encounter rates so that at
    # maximum size the species has feeding level f0
    if (is.null(crit_feed)) {
        crit_feed <- 0.7
    }
    f0 <- mizer::set_species_param_default(params@species_params, "f0", crit_feed)$f0

    # Get interaction of each species with detritus and algae
    ia <- params@species_params$interaction_algae
    id <- params@species_params$interaction_detritus

    # Calculate encounter rates divided by w^n of largest individuals
    E <- mizer::getEncounter(params)[, length(params@w)] /
        (params@w[length(params@w)]^n)

    # Calculate rho for each unstructured resource
    # f0*h/(1-f0) is the encounter rate when feeding level is f0
    # We subtract E so that if feeding level is too low, they eat
    # algae to replace it. For unstructured resources
    # encounter rate = rho * w^n * B_A, and multiply by interaction
    rho_alg <- pmax(0, f0 * params@species_params$h / (1 - f0) - E) * ia
    rho_det <- pmax(0, f0 * params@species_params$h / (1 - f0) - E) * id

    # # Set rho to 0 for predators with no max consumption rate
    # rho_alg[is.na(rho_alg)] <- 0
    # rho_det[is.na(rho_det)] <- 0

    # Store new rho values in species_params data frame
    params@species_params$rho_algae <- rho_alg
    params@species_params$rho_detritus <- rho_det

    # Calculate rho * w^n for use in algae and detritus dynamic functions
    rho_alg <- outer(params@species_params$rho_algae, params@w^0.86)
    rho_det <- outer(params@species_params$rho_detritus, params@w^n)

    # Store rho values in new params slots
    params@other_params$algae_params$rho_algae <- rho_alg
    params@other_params$detritus_params$rho_detritus <- rho_det

    ### Add components ----

    if (use_UR_cc == FALSE) {
        ### Algae Component - Add in algae ----
        params <- mizer::setComponent(
            params, "algae",
            initial_value = 1,
            dynamics_fun = "algae_dynamics",
            encounter_fun = "encounter_contribution",
            component_params = list(
                rho = rho_alg,
                capacity = params@other_params$algae_params$algae_capacity,
                growth = params@other_params$algae_params$algae_growth_initial
            )
        )

        ### Detritus component - Add in detritus ----
        params <- mizer::setComponent(
            params, "detritus",
            initial_value = 1,
            dynamics_fun = "detritus_dynamics",
            encounter_fun = "encounter_contribution",
            component_params = list(
                rho = rho_det,
                capacity = params@other_params$detritus_params$detritus_capacity,
                sen_decomp = params@other_params$detritus_params$sen_decomp,
                ext_decomp = params@other_params$detritus_params$ext_decomp,
                external = params@other_params$detritus_params$d_external
            )
        )
    } else if (use_UR_cc == TRUE) {
        ### Algae Component - Add in algae ----
        params <- mizer::setComponent(
            params, "algae",
            initial_value = 1,
            dynamics_fun = "algae_dynamics_cc",
            encounter_fun = "encounter_contribution",
            component_params = list(
                rho = rho_alg,
                capacity = params@other_params$algae_capacity,
                growth = params@other_params$initial_algae_growth
            )
        )

        ### Detritus component - Add in detritus ----
        params <- mizer::setComponent(
            params, "detritus",
            initial_value = 1,
            dynamics_fun = "detritus_dynamics_cc",
            encounter_fun = "encounter_contribution",
            component_params = list(
                rho = rho_det,
                sen_decomp = params@other_params$sen_decomp,
                ext_decomp = params@other_params$ext_decomp,
                capacity = params@other_params$detritus_capacity,
                external = params@other_params$external
            )
        )
    }

    # External mortality - Weight dependent ----
    if (include_ext_mort == TRUE) {
        ext_mort_params <- params@other_params[["ext_mort_params"]]

        # Change to allometric external mortality
        z0exp <- 1 - n
        nat_mort <- ext_mort_params$nat_mort
        nat_mort <- rep(nat_mort, nrow(species_params(params)))
        allo_mort <- outer(nat_mort, params@w^(z0exp))

        # Change the external mortality rate in the params object
        mizer::ext_mort(params) <- allo_mort
    } else {
        # Set coefficient for each species. Here we choose 0.1 for each species
        z0pre <- rep(z0pre, nrow(species_params(params)))
        z0exp <- 1 - n
        # Multiply by power of size with exponent, here chosen to be -1/4
        # The outer() function makes it an array species x size
        allo_mort <- outer(z0pre, w(params)^z0exp)

        # Change the external mortality rate in the params object
        mizer::ext_mort(params) <- allo_mort
    }

    # Register the extension chain and promote to the mizerReef S4 class ----
    # Rate overrides (Encounter, FeedingLevel, PredMort, Mort) are handled via
    # project*.mizerReef S3 methods defined in reef-project_methods.R, which
    # participate in the daisy-chain via NextMethod() rather than
    # setRateFunction(), making them composable with other extension packages.
    params@extensions <- mizer::getRegisteredExtensions()
    params <- mizer::coerceToExtensionClass(params)

    # Return object ----
    return(params)
}

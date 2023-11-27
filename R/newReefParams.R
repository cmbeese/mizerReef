#' Set up parameters for a mizerReef model
#'
#' Sets up a multi-species size spectrum model with additional unstructured
#' resource components, senescence mortality, and predation refuge.
#'
#' @inheritParams setURParams
#' @inheritParams setRho
#' @inheritParams setExtMortParams
#' @inheritParams setRefuge
#'
#' @param n Allometric growth exponent (also used as metabolic exponent p)
#' @param ... Extra parameters to be passed to [newMultispeciesParams()]
#'
#' @return An object of type \linkS4class{MizerParams}
#'
#' All the parameters will be mentioned in the following sections.
#' @inheritSection setURParams Adding unstructured resources
#' @inheritSection setExtMortParams Senescence mortality
#' @inheritSection setRefuge Setting the refuge profile
#'
#' @export
newReefParams <- function(species_params,
                          interaction = NULL,
                          min_w_pp = NA,
                          w_pp_cutoff = 1.0,
                          ext_mort_params = NULL,
                          UR_interaction,
                          algae_growth = NULL,
                          prop_decomp = NULL,
                          d.external = NULL,
                          crit_feed = NULL,
                          method,
                          min_ref_length = NULL,
                          max_protect = NULL,
                          tau = NULL,
                          method_params,
                          refuge_user = NULL,
                          bad_predator = NULL,
                          pisc = NULL,
                          n = 3/4, ...) {
    
    ## USE ORIGINAL MIZER FUNCTION TO INITIALIZE A MODEL
    params <- newMultispeciesParams(species_params = species_params,
                                    interaction = interaction,
                                    min_w_pp = min_w_pp,
                                    w_pp_cutoff = w_pp_cutoff,
                                    n = n, p = n, ...)
    
    # Check that a and b parameters are present for all species -
    # needed for l2w conversion
    if (anyNA(params@species_params[["a"]]) ||
        anyNA(params@species_params[["b"]])) {
        stop("There must be no NAs in the species_params columns 'a' and 'b'.")
    }
    
    ## ADD IN DETRITUS AND ALGAE
    params <- setURParams(params = params,
                          UR_interaction = UR_interaction,
                          algae_growth = algae_growth,
                          prop_decomp = prop_decomp,
                          d.external = d.external)
    
    # Determine the necessary detritus and algae encounter rates so that at
    #maximum size the group has feeding level f0
    if(is.null(crit_feed)){ crit_feed <- 0.6 }
    f0 <- set_species_param_default(params@species_params, "f0", crit_feed)$f0
    
    # Get interaction of each species with detritus and algae
    ia <- params@species_params$interaction_algae
    id <- params@species_params$interaction_detritus
    
    # Calculate encounter rates divided by w^n of largest individuals
    E <- getEncounter(params)[, length(params@w)] /
        (params@w[length(params@w)] ^ n)
    
    # Calculate rho for each unstructured resource
    # f0*h/(1-f0) is the encounter rate when feeding level is f0
    # We subtract E so that if feeding level is too low, they eat
    # algae to replace it. For unstructured resources
    # encounter rate = rho * w^n * B_A, and multiply by interaction
    rho_alg <- pmax(0, f0 * params@species_params$h / (1 - f0) - E) * ia
    rho_det <- pmax(0, f0 * params@species_params$h / (1 - f0) - E) * id
    
    # Store new rho values in species_params data frame
    params@species_params$rho_algae    <- rho_alg
    params@species_params$rho_detritus <- rho_det
    
    # # Function to calculate rho - why doesn't this work?
    # params <- setRho(params = params,
    #                  crit_feed = crit_feed, n = n)
 
    # Calculate rho * w^n for use in algae and detritus dynamic functions
    rho_alg <- outer(params@species_params$rho_algae, params@w ^ n)
    rho_det <- outer(params@species_params$rho_detritus, params@w ^ n)
    
    # Add in algae
    params <- setComponent(
        params, "algae", initial_value = 1,
        dynamics_fun = "algae_dynamics",
        encounter_fun = "encounter_contribution",
        component_params = list(rho = rho_alg,
                                algae_growth = params@other_params$algae_growth))
    
    # Add in detritus
    params <- setComponent(
        params, "detritus", initial_value = 1,
        dynamics_fun = "detritus_dynamics",
        encounter_fun = "encounter_contribution",
        component_params = list(rho = rho_det,
                                prop_decomp = params@other_params$prop_decomp,
                                d.external = params@other_params$d.external))

    ## ADD WEIGHT DEPENDENT MORTALITY
    params <- setExtMortParams(params = params,
                               ext_mort_params = ext_mort_params)

    ext_mort_params <- params@other_params[['ext_mort_params']]

    # Change to allometric external mortality
    z0exp     <- 1 - n
    nat_mort  <- ext_mort_params$nat_mort
    nat_mort  <- rep(nat_mort, nrow(species_params(params)))
    allo_mort <- outer(nat_mort, params@w^(z0exp))

    # Change the external mortality rate in the params object
    mizer::ext_mort(params) <- allo_mort

    # Add in senescence mortality
    params <- setRateFunction(params, "Mort", "reefMort")

    ## ADD REFUGE

    # Add parameters
    params <- setRefuge(params = params,
                        method = method,
                        min_ref_length = min_ref_length,
                        max_protect = max_protect,
                        tau = tau,
                        method_params = method_params,
                        refuge_user = refuge_user,
                        bad_predator = bad_predator,
                        pisc = pisc, ...)

    # Replace mizerRate functions with mizerReef versions
    # mizerRates
    params <- setRateFunction(params, "Rates", "reefRates")
    # mizerEncounter
    params <- setRateFunction(params, "Encounter", "reefEncounter")
    # mizerPredMort
    params <- setRateFunction(params, "PredMort", "reefPredMort")

    params
}

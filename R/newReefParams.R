#' Set up parameters for a mizerReef model
#'
#' Sets up a multi-species size spectrum model with additional unstructured
#' resource components, senescence mortality, and predation refuge.
#' 
#' @inheritSection setURParams Adding unstructured resources
#' @inheritSection setExtMortParams Senescence mortality
#' @inheritSection setRefuge Setting the refuge profile
#' 
#' @export
#' @inheritParams setURParams
#' @inheritParams setExtMortParams
#' @inheritParams setRefuge
#' @param species_params A functional group parameter data frame
#' @param interaction The group specific interaction matrix, \eqn{\theta_{ij}}
#' @param min_w_pp Minimum size of plankton in grams
#' @param w_pp_cutoff Maximum size of plankton in grams
#' @param n Growth exponent (also used as metabolic exponent p)
#' @param crit_feed Critical feeding level
#' @param ... Extra parameters to be passed to [newMultispeciesParams()]
#'
#' @param n Allometric growth exponent (also used as metabolic exponent p)
#' @param ... Extra parameters to be passed to [newMultispeciesParams()]
#' @concept setup
#' @return An object of type \linkS4class{MizerParams}
newReefParams <- function(# Original mizer parameters
                            species_params, interaction = NULL, 
                            crit_feed = NULL,
                            min_w_pp = NA, w_pp_cutoff = 1.0, n = 3/4,
                          # Parameters for setting up refuge
                            method, method_params, 
                            refuge_user = NULL, bad_pred = NULL, 
                            satiation = NULL,
                            a_bar = NULL, b_bar = NULL,
                            w_settle = NULL, max_protect = NULL, tau = NULL,
                          # Parameters for unstructured resources
                            UR_interaction, exp_alg = NULL, exp_det = NULL,
                            algae_growth = NULL, 
                            scale_rho_a = NULL, scale_rho_d = NULL,
                            prop_decomp = NULL, d.external = NULL,
                          # Parameters for external mortality
                            ext_mort_params = NULL, ...) {
    
    ## Initialize model with newMultispeciesParams ----
    params <- newMultispeciesParams(species_params = species_params,
                                    interaction = interaction,
                                    min_w_pp = min_w_pp,
                                    w_pp_cutoff = w_pp_cutoff,
                                    n = n, p = n, ...)
    
    # Add parameters ----
    ### Refuge ----
    params <- setRefuge(params = params, method = method, 
                        method_params = method_params,
                        refuge_user = refuge_user, 
                        bad_pred = bad_pred,
                        satiation = satiation,
                        a_bar = a_bar, b_bar = b_bar,
                        w_settle = w_settle, 
                        max_protect = max_protect, tau = tau, ...)
    
    params <- getRefuge(params,...)
    
    if(is.null(exp_det)){exp_det <- n}
    ### Unstructured resources ----
    params <- setURParams(params = params,
                          UR_interaction = UR_interaction,
                          exp_alg = exp_alg,
                          exp_det = exp_det,
                          algae_growth = algae_growth,
                          scale_rho_a = scale_rho_a,
                          scale_rho_d = scale_rho_d,
                          prop_decomp = prop_decomp,
                          d.external = d.external)
    
    ### External mortality ----
    params <- setExtMortParams(params = params,
                               ext_mort_params = ext_mort_params)
    
    
    # Algae & Detritus ----
    ### Calculate Rho ----
        # Determine the necessary detritus and algae encounter rates so that at
        # maximum size the group has feeding level f0
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
        
        # Set rho to 0 for predators with no max consumption rate
        rho_alg[is.na(rho_alg)] <- 0
        rho_det[is.na(rho_det)] <- 0
        
        # Pull scaling values
        scale_rho_a <- params@other_params$scale_rho_a
        scale_rho_d <- params@other_params$scale_rho_d
        
        # Store new rho values in species_params data frame
        params@species_params$rho_algae    <- scale_rho_a*rho_alg
        params@species_params$rho_detritus <- scale_rho_d*rho_det
    
        # Calculate rho * w^n for use in algae and detritus dynamic functions
        exp_alg <- params@other_params$exp_alg
        exp_det <- params@other_params$exp_det
        rho_alg <- outer(params@species_params$rho_algae, params@w ^ exp_alg)
        rho_det <- outer(params@species_params$rho_detritus, params@w ^ exp_det)
    
    ### Algae Component - Add in algae ----
    params <- setComponent(
        params, "algae", initial_value = 1,
        dynamics_fun = "algae_dynamics",
        encounter_fun = "encounter_contribution",
        component_params = list(rho = rho_alg,
                                algae_growth = params@other_params$algae_growth))
    
    ### Detritus component - Add in detritus ----
    params <- setComponent(
        params, "detritus", initial_value = 1,
        dynamics_fun = "detritus_dynamics",
        encounter_fun = "encounter_contribution",
        component_params = list(rho = rho_det,
                                prop_decomp = params@other_params$prop_decomp,
                                d.external  = params@other_params$d.external))

    # External mortality - Weight dependent ----
        ext_mort_params <- params@other_params[['ext_mort_params']]
    
        # Change to allometric external mortality
        z0exp     <- 1 - n
        nat_mort  <- ext_mort_params$nat_mort
        nat_mort  <- rep(nat_mort, nrow(species_params(params)))
        allo_mort <- outer(nat_mort, params@w^(z0exp))
    
        # Change the external mortality rate in the params object
        mizer::ext_mort(params) <- allo_mort

    # Changes rates for refuge ----
        # Replace mizerRate functions with mizerReef versions
        # mizerRates
        params <- setRateFunction(params, "Rates", "reefRates")
        # mizerEncounter
        params <- setRateFunction(params, "Encounter", "reefEncounter")
        # mizerFeedingLevel
        params <- setRateFunction(params, "FeedingLevel", "reefFeedingLevel")
        # mizerPredMort
        params <- setRateFunction(params, "PredMort", "reefPredMort")
        # Add in senescence mortality
        params <- setRateFunction(params, "Mort", "reefMort")
        
        # Change scale model function
        params <- customFunction("scaleModel", reefScaleModel)
    
    # Return object ----    
    return(params)
}

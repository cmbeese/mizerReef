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
#' @param group_params  A functional group parameter data frame containing at
#'                      least the name of each functional group, their
#'                      observed abundances, and the cut-off size for 
#'                      observations in grams.
#' @param interaction The group specific interaction matrix, \eqn{\theta_{ij}}
#' 
#' @param carry_capacity A boolean value that indicates whether the user wants
#'                      to implement a carrying capacity for unstructured 
#'                      resources. Defaults to FALSE.
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
                            group_params, interaction = NULL, 
                            crit_feed = NULL,
                            min_w_pp = NA, w_pp_cutoff = 1, 
                            n = 0.75,
                          # Parameters for setting up refuge
                            method, method_params, 
                            refuge_user = NULL, bad_pred = NULL, 
                            satiation = NULL,
                            a_bar = NULL, b_bar = NULL,
                            w_settle = NULL, max_protect = NULL, 
                            tau = NULL,
                          # Parameters for unstructured resources
                            UR_interaction,
                            carry_capacity = FALSE,
                            initial_algae_growth = NULL, 
                            algae_capacity = NULL,
                            detritus_capacity = NULL,
                            sen_decomp = NULL, ext_decomp = NULL, 
                            initial_d_external = NULL,
                          # Parameters for external mortality
                            ext_mort_params = NULL,
                            include_ext_mort = TRUE,
                            include_sen_mort = TRUE,
                            z0pre = 0.2, ...) {
    
    ## Initialize model with newMultispeciesParams ----
    params <- newMultispeciesParams(species_params = group_params,
                                    interaction = interaction,
                                    min_w_pp = min_w_pp,
                                    w_pp_cutoff = w_pp_cutoff,
                                    n = n, p = n, ...)
    
    # Change resource color
    params@linecolour["Resource"] <- "lightseagreen"
    params@linecolour["Fishing"] <- "blue"
    
    # Save include_sen_mort
    params@other_params$include_sen_mort <- include_sen_mort
    
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
    
    ### Unstructured resources ----
    params <- setURParams(params = params,
                          UR_interaction = UR_interaction,
                          initial_algae_growth = initial_algae_growth,
                          carry_capacity = carry_capacity,
                          algae_capacity = algae_capacity,
                          detritus_capacity = detritus_capacity,
                          sen_decomp = sen_decomp, 
                          ext_decomp = ext_decomp, 
                          initial_d_external = initial_d_external)
    
    ### External mortality ----
    params <- setExtMortParams(params = params,
                               ext_mort_params = ext_mort_params)
    
    
    # Algae & Detritus ----
    ### Calculate Rho ----
        # Determine the necessary detritus and algae encounter rates so that at
        # maximum size the group has feeding level f0
        if(is.null(crit_feed)){ crit_feed <- 0.7 }
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
        
        # # Set rho to 0 for predators with no max consumption rate
        # rho_alg[is.na(rho_alg)] <- 0
        # rho_det[is.na(rho_det)] <- 0
        
        # Store new rho values in species_params data frame
        params@species_params$rho_algae    <- rho_alg
        params@species_params$rho_detritus <- rho_det
    
        # Calculate rho * w^n for use in algae and detritus dynamic functions
        rho_alg <- outer(params@species_params$rho_algae, params@w ^ 0.86)
        rho_det <- outer(params@species_params$rho_detritus, params@w ^ n)
    
    if (carry_capacity == FALSE){
        ### Algae Component - Add in algae ----
        params <- setComponent(
            params, "algae", initial_value = 1,
            dynamics_fun = "algae_dynamics",
            encounter_fun = "encounter_contribution",
            component_params = list(rho = rho_alg,
                                    growth = params@other_params$initial_algae_growth))
        
        ### Detritus component - Add in detritus ----
        params <- setComponent(
            params, "detritus", initial_value = 1,
            dynamics_fun = "detritus_dynamics",
            encounter_fun = "encounter_contribution",
            component_params = list(rho = rho_det,
                                    sen_decomp = params@other_params$sen_decomp,
                                    ext_decomp = params@other_params$ext_decomp,
                                    external   = params@other_params$initial_d_external))
        
    } else if (carry_capacity == TRUE){
        ### Algae Component - Add in algae ----
        params <- setComponent(
            params, "algae", initial_value = 1,
            dynamics_fun = "algae_dynamics_cc",
            encounter_fun = "encounter_contribution",
            component_params = list(rho = rho_alg,
                                    capacity = params@other_params$algae_capacity,
                                    growth   = params@other_params$initial_algae_growth))
        
        ### Detritus component - Add in detritus ----
        params <- setComponent(
            params, "detritus", initial_value = 1,
            dynamics_fun = "detritus_dynamics_cc",
            encounter_fun = "encounter_contribution",
            component_params = list(rho = rho_det,
                                    sen_decomp = params@other_params$sen_decomp,
                                    ext_decomp = params@other_params$ext_decomp,
                                    capacity   = params@other_params$detritus_capacity,
                                    external   = params@other_params$initial_d_external))
    }

    # External mortality - Weight dependent ----
    if (include_ext_mort == TRUE){
        ext_mort_params <- params@other_params[['ext_mort_params']]
    
        # Change to allometric external mortality
        z0exp     <- 1 - n
        nat_mort  <- ext_mort_params$nat_mort
        nat_mort  <- rep(nat_mort, nrow(species_params(params)))
        allo_mort <- outer(nat_mort, params@w^(z0exp))
    
        # Change the external mortality rate in the params object
        mizer::ext_mort(params) <- allo_mort
    } else {
        # Set coefficient for each species. Here we choose 0.1 for each species
        z0pre <- rep(z0pre, nrow(species_params(params)))
        z0exp     <- 1 - n
        # Multiply by power of size with exponent, here chosen to be -1/4
        # The outer() function makes it an array species x size
        allo_mort <- outer(z0pre, w(params)^z0exp)
        
        # Change the external mortality rate in the params object
        mizer::ext_mort(params) <- allo_mort
    }

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
    
    # Return object ----    
    return(params)
}

#' Calculate fisheries productivity for each functional group above a chosen
#' minimum fishing size
#'
#' @section Potential fisheries productivity:
#' 
#' Fisheries productivity refers to the rate at which fish biomass is produced
#' and available for harvest in a given area over a given period of time.
#' Productivity cannot be measured in situ.
#'
#' The productivity \eqn{P_i(w)} of functional group \eqn{i} is given by
#' 
#     \deqn{P_i(w) = \int_w^{w+dw} \left( N_i(w) + g_i(w) \right) w \, dw.}
#          {P_i(w) = \int_w^{w+dw} \left( N_i(w) + g_i(w) \right) w \, dw.}
#'          
#'     \deqn{P_i(w) = \int_w^{w+dw} \left( N_i(w) + E_{R.i}(w) \right) w \, dw.}
#'          {P_i(w) = \int_w^{w+dw} ( N_i(w) + E_{R.i}(w)) w \, dw.}
#'
#'  \eqn{N_i(w)} is the abundance density (1/m^-2) and \eqn{E_{R.i}} 
#  \eqn{g_i(w)} 
#'  is the energy rate available for growth after metabolism and movement have 
#'  been accounted for (grams/year).
#'
#' The productivity is calculated for all fish in the size range between
#' `min_fishing_length` and `max_fishing_length`which can be the same for 
#' all functional groups or can be specified as a vector with one value for 
#' each species in the model.`min_fishing_length` defaults to \eqn{7 cm} 
#' regardless of functional group and `max_fishing_length` defaults to the 
#' maximum weight in the model. 
#'
#' @param object    An object of class `MizerParams` or `MizerSim`. If given a
#'                  \linkS4class{MizerSim} object, uses the growth rates at the 
#'                  final time of a simulation to calculate productivity. If 
#'                  given a \linkS4class{MizerParams} object, uses the initial 
#'                  growth rates.
#'                  
#@param total A boolean value that indicates whether the user wants the
#             mass specific productivity for each functional group (species
#             by size) or the total productivity for each functional group.
#             Default value is true.
#'              
#' @param min_fishing_l The minimum length (cm) of fished individuals for
#'                      productivity estimates. Defaults to 7 cm.
#'                      
#' @param max_fishing_l The maximum length (cm) of fished individuals for
#'                      productivity estimates. Defaults to max length.
#'                       
#' @inheritDotParams mizer::get_size_range_array -params
#'
#' @return If called with a MizerParams object, a vector with the productivity
#'   in grams/year/m^-2 for each functional group in the model. If called with
#'   a MizerSim object, an array (time x species) containing the productivity
#'   at each time step for all species.
#'
#' @export
#' @family summary functions
#' @concept summary
getProductivity <- function(object,
                            min_fishing_l = NULL, 
                            max_fishing_l = NULL,...) {
    
    if (is(object, "MizerSim")) {
        stop('This functionality is not set up yet.')
    }

    if (is(object, "MizerParams")) {
        params <- object
        
        if(is.null(min_fishing_l)){ min_fishing_l <- 7 }
        if(is.null(max_fishing_l)){ max_fishing_l <- max(params@species_params$l_max) }
        
        # why isnt this working?
        size_range <- mizer::get_size_range_array(params, 
                                                  min_l = min_fishing_l, 
                                                  max_l = max_fishing_l, ...)

        # pr <- mizer::getEGrowth(params)
        pr <- mizer::getEReproAndGrowth(params)
        prod <- (pr * params@initial_n * size_range) %*% (params@w * params@dw)
        # This seems like growth times biomass?
        # Alice's old code:
        # prod_pred <- colSums(10^x[ref:end]*GG.u[ref:end,]*U[ref:end,]*dx)
        # 10^x[ref:end] is weights of all fish over 7 cm in length
        # growth of predators over 7 cm in length GG.u[ref:end,]
        # abundance of predators over 7 cm in length U[ref:end,]
        # dx is width of the size bins
        # So for each species prod <- colSums(w*g*n*dw)
        
        # if(total == FALSE) {return(pr)}
        # if(total == TRUE){
        # Not sure which one it's supposed to be
        # prod <- (pr * params@w) %*% params@dw
        # prod <- pr  %*% (params@w * params@dw)
        return(prod[, , drop = TRUE])
        #}
    }
    stop("'object' should be a MizerParams or a MizerSim object")
}




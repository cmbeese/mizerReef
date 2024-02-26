#' Calculate fisheries productivity for each species group
#' 
#' This function calculates the total potential fisheries productivity by
#' species group in a given size range using species abundances and
#' the [getEGrowthTime()] function.
#'
#' @section Potential fisheries productivity:
#' 
#'  Productivity refers to the rate at which fish biomass is produced
#'  and available for harvest in a given area over a given period of time.
#'  Productivity cannot be measured in situ.
#'      
#'  The productivity \eqn{P_i(w)} of species group \eqn{i} is given by
#'  
#'  \deqn{P_i(w) = \int_w^{w+dw} \left( N_i(w) + g_i(w) \right) w \, dw.}{
#'       P_i(w) = \int_w^{w+dw} ( N_i(w) + g_i(w) ) w \, dw.}
#'      
#'  \eqn{N_i(w)} is the abundance density \eqn{(no./m^{2})} and 
#'  \eqn{g_i(w)} is the energy rate available for growth after metabolism 
#'  and movement have been accounted for \eqn{(grams/year)}.
#
#'  The productivity is calculated for all fish in the size range between
#'  `min_fishing_length` and `max_fishing_length`. These lengths can be the 
#'  same for all groups or can be specified as a vector with one value for 
#'  each species in the model. The minimum length defaults to \eqn{7 cm} 
#'  regardless of species group and maximum length defaults to the 
#'  maximum weight in the model. 
#'
#' @param object    An object of class `MizerParams` or `MizerSim`. If given a
#'                  \linkS4class{MizerSim} object, uses the growth rates at the 
#'                  final time of a simulation to calculate productivity. If 
#'                  given a \linkS4class{MizerParams} object, uses the initial 
#'                  growth rates.
#'                  
#' @param include_repro A boolean value that indicates whether to include
#'                      energy for reproduction in productivity estimates.
#'                      Defaults to using [getEGrowthTime()] if FALSE or
#'                      mizer's [mizer::getEReproAndGrowth()] if TRUE.
#'              
#' @param min_fishing_l The minimum length (cm) of fished individuals for
#'                      productivity estimates. Defaults to 7 cm.
#'                      
#' @param max_fishing_l The maximum length (cm) of fished individuals for
#'                      productivity estimates. Defaults to max length.
#'                      
#' @param time_range The time range (either a vector of values, a vector of min
#'                   and max time, or a single value) to provide productivity
#'                   for. Default is the final time step. Ignored when called 
#'                   with a \linkS4class{MizerParams} object.
#'                   
#' @inheritDotParams mizer::get_size_range_array
#'
#' @return  If called with a MizerParams object, a vector with the productivity
#'          in \eqn{\frac{grams}{m^2}/year} for each species group in the model. 
#'          If called with a MizerSim object, an array (time x species) 
#'          containing the total productivity at each time step for each 
#'          species.
#' 
#' @importFrom plyr aaply
#' @export
#' @family summary functions
#' @concept summary
#' @seealso [getEGrowthTime()],[plotProductivity()]
getProductivity <- function(object,
                            time_range = NULL,
                            include_repro = FALSE,
                            min_fishing_l = NULL, 
                            max_fishing_l = NULL,...) {
    
    if (is(object, "MizerParams")) {
        
        params <- object
        assert_that(is(params, "MizerParams"))
        
        # Set default fishing sizes if not provided
        if(is.null(min_fishing_l)){ min_fishing_l <- 7 }
        if(is.null(max_fishing_l)){ 
            max_fishing_l <- max(params@species_params$l_max) 
        }
        
        # Get matrix of true false values for fish in size range
        size_range <- mizer::get_size_range_array(params, 
                                                  min_l = min_fishing_l, 
                                                  max_l = max_fishing_l, ...)
        
        if (include_repro == FALSE){ energy <- getEGrowthTime(params) }
        if (include_repro == TRUE) { 
            energy <- mizer::getEReproAndGrowth(params) 
        }
        
        prod <- ((energy * params@initial_n * size_range) %*% 
                     params@dw)[,,drop = TRUE]
        
        return(prod)
    
    } else if (is(object, "MizerSim")) {
        
        # If no time range is given, default to the final time step
        sim <- object
        if (missing(time_range)) {
            time_range <- max(as.numeric(dimnames(sim@n)$time))
        }
        
        # Get matrix of true false values for fish in size range
        size_range <- mizer::get_size_range_array(sim@params, 
                                                  min_l = min_fishing_l, 
                                                  max_l = max_fishing_l,...)
        
        # Get matrix of true false values for times in time range
        time_elements <- mizer::get_time_elements(sim, time_range)
        
        if (include_repro == FALSE){ e_time <- getEGrowthTime(sim, time_range) }
        if (include_repro == TRUE) { 
            stop('This functionality is not set up yet.')
        }
        
        prod <- plyr::aaply(which(time_elements), 1, function(x) {
            # Necessary as we only want single time step but may only have 1
            # species which makes using drop impossible
            n <- array(sim@n[x, , ], dim = dim(sim@n)[2:3])
            dimnames(n) <- dimnames(sim@n)[2:3]
            energy <- array(e_time[,,x], dim = dim(e_time)[2:3])
            dimnames(energy) <- dimnames(e_time)[2:3]
            prod <- ((energy * n * size_range) %*% params@dw)[,, drop = TRUE]
            return(prod)
            }, .drop = FALSE)
            
        return(prod)
    } else {
        stop("Object should be a MizerParams or MizerSim object.")
    }
}




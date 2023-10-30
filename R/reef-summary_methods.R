#' Calculate the total fisheries productivity of each functional group within
#' a size range at each time step
#'
#' Fisheries productivity refers to the rate at which fish biomass is produced
#' and available for harvest in a given area over a given period of time.
#' Productivity cannot be measured in situ.
#'
#' The productivity \eqn{P_i(w)} of functional group \eqn{i} is given by
#'
#'  \deqn{P_i(w) = \int_w^{w+dw} N_i(w) g_i(w) dw}
#'          {P_i(w) = int_w^{w+dw} N_i(w) g_i(w) dw}
#'
#'  \eqn{N_i(w)} is the abundance density (1/m^-2) and \eqn{g_i(w)} is the
#'  energy rate available for growth after metabolism, movement and
#'  reproduction have been accounted for (grams/year).
#'
#' The productivity is calculated for all fish larger than input parameter
#' `min_fishing_length`which defaults to \eqn{7 cm} regardless of functional
#' group.
#'
#' @param object An object of class `MizerParams` or `MizerSim`.If given a
#'  \linkS4class{MizerSim} object, uses the growth rates at the final time of a
#'   simulation to calculate productivity. If given a \linkS4class{MizerParams}
#'   object, uses the initial growth rates.
#' @param min_fishing_l The minimum size of fished inidividuals for
#'      productivity estimates. Defaults to 7 cm.
#'
#' @return If called with a MizerParams object, a vector with the productivity
#'   in grams/year/m^-2 for each functional group in the model. If called with
#'   a MizerSim object, an array (time x species) containing the productivity
#'   at each time step for all species.
#'
#' @export
#' @family summary functions
#' @concept summary_function
getProductivity <- function(object,
                            min_fishing_l = 7,...) {
    if (is(object, "MizerSim")) {
        stop('This functionality is not set up yet you dumbass.')
    }

    if (is(object, "MizerParams")) {
        params <- object
        size_range <- get_size_range_array(params,
                                           min_l = min_fishing_l,...)

        # NEED TO CHECK ON HOW THIS IS CALCULATED!
        g <- mizer::getEGrowth(params)
        prod <- params@initial_n * size_range * g
        return((prod %*% (params@w * params@dw))[, , drop = TRUE])
    }

    stop("'object' should be a MizerParams or a MizerSim object")
}

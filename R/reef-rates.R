#' Get vulnerability level at in time range t
#'
#' Returns the proportion of fish at size \eqn{w} that are not hidden in
#' predation refuge and thus vulnerable to being encountered by predators.
#' 
#' This function uses [reefVulnerable()] to calculate the vulnerability to 
#' predation.
#' 
#' @inherit reefVulnerable
#' 
#' @param object A `MizerParams` object or a `MizerSim` object
#' 
#' @inheritParams reefRates
#' 
#' @inheritParams mizer::get_time_elements
#' 
#' @param drop  If `TRUE` then any dimension of length 1 will be removed
#'              from the returned array.
#'
#' @return  If a `MizerParams` object is passed in, the function returns a two
#'          dimensional array (prey species x prey size) based on the
#'          abundances also passed in.
#'   
#'          If a `MizerSim` object is passed in, the function returns a three
#'          dimensional array (time step x prey species x prey size) 
#'          with the vulnerability calculated at every time step in the 
#'          simulation. If \code{drop = TRUE} then the dimension of length 1 
#'          will be removed from the returned array.
#' 
#' @import plyr 
#' @export
#' @concept refugeRates
#' @family rate functions
getVulnerable <- function(object, n, n_pp, n_other,
                          time_range, drop = TRUE, ...) {
    
    if (is(object, "MizerParams")) {
        # if params -----
        params <- mizer::validParams(object)
        if (missing(time_range)) time_range <- 0
        t <- min(time_range)
        if (missing(n)) n <- params@initial_n
        if (missing(n_pp)) n_pp <- params@initial_n_pp
        if (missing(n_other)) n_other <- params@initial_n_other
        # calculate vulnerability
        vulnerable <- reefVulnerable(params, n = n, n_pp = n_pp, 
                                     n_other = n_other, t = t)
        dimnames(vulnerable) <- dimnames(params@metab)
        return(vulnerable)
        
    } else {
        # if object ----
        sim <- object
        if (missing(time_range)) {
            time_range <- dimnames(sim@n)$time
        }
        time_elements <- mizer::get_time_elements(sim, time_range)
        vul_time <- plyr::aaply(which(time_elements), 1, function(x) {
            # Necessary as we only want single time step but may only have 1
            # species which makes using drop impossible
            n <- array(sim@n[x, , ], dim = dim(sim@n)[2:3])
            dimnames(n) <- dimnames(sim@n)[2:3]
            n_other <- sim@n_other[x, ]
            names(n_other) <- dimnames(sim@n_other)$component
            t <- as.numeric(dimnames(sim@n)$time[[x]])
            vul <- getVulnerable(sim@params, n = n,
                                 n_pp = sim@n_pp[x, ],
                                 n_other = n_other,
                                 time_range = t)
            return(vul)
        }, .drop = FALSE)
        # Before we drop dimensions we want to set the time dimname
        names(dimnames(vul_time))[[1]] <- "time"
        vulnerable <- vul_time[, , , drop = drop]
        return(vulnerable)
    }
}



#' Get the size specific senescence mortality rate
#'
#' Returns the rate of senescence mortality at each size by functional group.
#'
#' @inherit reefSenMort
#'
#' @export
#' @concept extmort
#' @family rate functions
getSenMort <- function(params, n = initialN(params),
                       n_pp = initialNResource(params),
                       n_other = initialNOther(params),
                       t = 0, ...) {
    
    params <- validParams(params)
    assert_that(is.array(n),
                is.numeric(n_pp),
                is.list(n_other),
                splus2R::is.number(t),
                identical(dim(n), dim(params@initial_n)),
                identical(length(n_pp), length(params@initial_n_pp)),
                identical(length(n_other), length(params@initial_n_other))
    )
    
    sen_mort <- reefSenMort(params, n = n, n_pp = n_pp, 
                            n_other = n_other, t = t)
    sen_mort
}


#' Get energy rate available for growth through time
#'
#' Calculates the energy rate \eqn{g_i(w)} (grams/year) available by 
#' species and size for growth after metabolism, movement and 
#' reproduction have been accounted for.
#' 
#' @inheritParams reefRates
#'   
#' @return If a `MizerParams` object is passed in, the function returns a two
#'   dimensional array (predator species x predator size) based on the
#'   abundances also passed in.
#'   If a `MizerSim` object is passed in, the function returns a three
#'   dimensional array (time step x predator species x predator size) with the
#'   energy for growth calculated at every time step in the simulation.
#'   If \code{drop = TRUE} then the dimension of length 1 will be removed from
#'   the returned array.
#' @export
#' @concept summary
getEGrowthTime <- function(params, n, n_pp, n_other,
                           time_range,
                           drop = FALSE, ...) {
    
    if (is(object, "MizerParams")) {
        params <- validParams(params)
        f <- get(params@rates_funcs$EGrowth)
        
        # Get any missing arguments
        if (missing(time_range)) time_range <- 0
        t <- min(time_range)
        if (missing(n)) n <- params@initial_n
        if (missing(n_pp)) n_pp <- params@initial_n_pp
        if (missing(n_other)) n_other <- params@initial_n_other
        
        # Calculate growth
        g <- f(params, n = n, n_pp = n_pp, n_other = n_other, t = t, 
               e_repro = getERepro(params, n = n, n_pp = n_pp, 
                                   n_other = n_other, t = t), 
               e = getEReproAndGrowth(params, n = n, n_pp = n_pp, 
                                      n_other = n_other, t = t))
        dimnames(g) <- dimnames(params@metab)
        
        return(g)
        
    } else { 

        sim <- object
        if (missing(time_range)) {
            time_range <- dimnames(sim@n)$time
        }
        time_elements <- mizer::get_time_elements(sim, time_range)
        grow_time <- plyr::aaply(which(time_elements), 1, function(x) {
            # Necessary as we only want single time step but may only have 1
            # species which makes using drop impossible
            n <- array(sim@n[x, , ], dim = dim(sim@n)[2:3])
            dimnames(n) <- dimnames(sim@n)[2:3]
            n_other <- sim@n_other[x, ]
            names(n_other) <- dimnames(sim@n_other)$component
            t <- as.numeric(dimnames(sim@n)$time[[x]])
            grow <- getEGrowthTime(sim@params, n = n,
                                   n_pp = sim@n_pp[x, ],
                                   n_other = n_other,
                                   time_range = t)
            return(grow)
        }, .drop = FALSE)
    
    # Before we drop dimensions we want to set the time dimname
    names(dimnames(grow_time))[[1]] <- "time"
    grow_time <- grow_time[, , , drop = drop]
    return(grow_time)
    }
}
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
#' @return  If a `MizerParams` object is passed in, the function returns a
#'          numeric vector of the refuge density for each size bin.
#'
#'          If a `MizerSim` object is passed in, the function returns a two
#'          dimensional array (time step x refuge size bin) with the refuge
#'          density calculated at every time step in the simulation. If
#'          \code{drop = TRUE} then the dimension of length 1 will be
#'          removed from the returned array.
#'
#' @export
#' @concept degradation
#' @family rate functions
getDegrade <- function(object, n, n_pp, n_other,
                       time_range, drop = TRUE, ...) {
    if (inherits(object, "MizerParams")) {
        # params -----
        params <- mizer::validParams(object)
        if (missing(time_range)) time_range <- 0
        t <- min(time_range)
        if (missing(n)) n <- params@initial_n
        if (missing(n_pp)) n_pp <- params@initial_n_pp
        if (missing(n_other)) n_other <- params@initial_n_other

        # calculate vulnerability
        degrade <- reefDegrade(params,
            n = n, n_pp = n_pp,
            n_other = n_other, t = t
        )

        return(degrade)
    } else {
        # sim ----
        sim <- object
        if (missing(time_range)) {
            time_range <- dimnames(sim@n)$time
        }
        time_elements <- mizer::get_time_elements(sim, time_range)
        deg_time <- plyr::aaply(which(time_elements), 1, function(x) {
            # Necessary as we only want single time step but may only have 1
            # species which makes using drop impossible
            n <- array(sim@n[x, , ], dim = dim(sim@n)[2:3])
            dimnames(n) <- dimnames(sim@n)[2:3]
            n_other <- sim@n_other[x, ]
            names(n_other) <- dimnames(sim@n_other)$component
            t <- as.numeric(dimnames(sim@n)$time[[x]])
            deg <- getDegrade(sim@params,
                n = n, n_pp = sim@n_pp[x, ],
                n_other = n_other, time_range = t
            )
            return(deg)
        }, .drop = FALSE)
        # Before we drop dimensions we want to set the time dimname
        names(dimnames(deg_time))[[1]] <- "time"
        # reefDegrade() returns a 1D vector (refuge density by size bin), so
        # aaply()-ing it over time steps gives a 2D [time, size bin] array,
        # not 3D - indexing with a third comma here always errored before
        # this was caught (no test previously exercised this MizerSim path).
        degrade <- deg_time[, , drop = drop]
        return(degrade)
    }
}


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
#' @export
#' @concept refugeRates
#' @family rate functions
getVulnerable <- function(object, n, n_pp, n_other,
                          time_range, drop = TRUE, ...) {
    if (inherits(object, "MizerParams")) {
        # params -----
        params <- mizer::validParams(object)
        if (missing(time_range)) time_range <- 0
        t <- min(time_range)
        if (missing(n)) n <- params@initial_n
        if (missing(n_pp)) n_pp <- params@initial_n_pp
        if (missing(n_other)) n_other <- params@initial_n_other

        new_rd <- getDegrade(params,
            n = n, n_pp = n_pp,
            n_other = n_other, time_range = t
        )

        # calculate vulnerability
        vulnerable <- reefVulnerable(params,
            n = n, n_pp = n_pp,
            n_other = n_other, t = t,
            new_rd = new_rd
        )
        dimnames(vulnerable) <- dimnames(params@metab)
        return(vulnerable)
    } else {
        # sim ----
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
            new_rd <- getDegrade(sim@params,
                n = n,
                n_pp = sim@n_pp[x, ],
                n_other = n_other, time_range = t
            )
            vul <- getVulnerable(sim@params,
                n = n,
                n_pp = sim@n_pp[x, ],
                n_other = n_other,
                time_range = t,
                new_rd = new_rd
            )
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
#' @inheritParams reefRates
#'
#' @export
#' @concept extmort
#' @family rate functions
getSenMort <- function(params, n = initialN(params),
                       n_pp = params@initial_n_pp,
                       n_other = initialNOther(params),
                       t = 0, ...) {
    params <- validParams(params)
    assert_that(
        is.array(n),
        is.numeric(n_pp),
        is.list(n_other),
        is.number(t),
        identical(dim(n), dim(params@initial_n)),
        identical(length(n_pp), length(params@initial_n_pp)),
        identical(length(n_other), length(params@initial_n_other))
    )

    sen_mort <- reefSenMort(params,
        n = n, n_pp = n_pp,
        n_other = n_other, t = t
    )
    sen_mort
}


#' Get energy rate available for growth through time
#'
#' Calculates the energy rate \eqn{g_i(w)} (grams/year) available by
#' species and size for growth after metabolism, movement and
#' reproduction have been accounted for.
#'
#' @param object A `MizerParams` object or a `MizerSim` object
#'
#' @param drop If \code{drop = TRUE} then the dimension of length 1 will be
#'      removed from the returned array.
#'
#' @inheritParams reefRates
#'
#' @inheritParams mizer::get_time_elements
#'
#' @return If a `MizerParams` object is passed in, the function returns a two
#'   dimensional array (predator species x predator size) based on the
#'   abundances also passed in.
#'   If a `MizerSim` object is passed in, the function returns a three
#'   dimensional array (time step x predator species x predator size) with the
#'   energy for growth calculated at every time step in the simulation.
#'   If \code{drop = TRUE} then the dimension of length 1 will be removed from
#'   the returned array.
#'
#' @export
#' @concept summary
#' @seealso [getProductivity()]
getEGrowthTime <- function(object, n, n_pp, n_other,
                           time_range,
                           drop = FALSE, ...) {
    if (inherits(object, "MizerParams")) {
        params <- object
        params <- validParams(params)
        f <- get(params@rates_funcs$EGrowth)

        # Get any missing arguments
        if (missing(time_range)) time_range <- 0
        t <- min(time_range)
        if (missing(n)) n <- params@initial_n
        if (missing(n_pp)) n_pp <- params@initial_n_pp
        if (missing(n_other)) n_other <- params@initial_n_other

        # Calculate growth
        g <- f(params,
            n = n, n_pp = n_pp, n_other = n_other, t = t,
            e_repro = getERepro(params,
                n = n, n_pp = n_pp,
                n_other = n_other, t = t
            ),
            e = getEReproAndGrowth(params,
                n = n, n_pp = n_pp,
                n_other = n_other, t = t
            )
        )
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
            grow <- getEGrowthTime(sim@params,
                n = n,
                n_pp = sim@n_pp[x, ],
                n_other = n_other,
                time_range = t
            )
            return(grow)
        }, .drop = FALSE)

        # Before we drop dimensions we want to set the time dimname
        names(dimnames(grow_time))[[1]] <- "time"
        grow_time <- grow_time[, , , drop = drop]
        return(grow_time)
    }
}


#' Get the diet composition of a mizerReef model
#'
#' Extends [mizer::getDiet()] so that the diet agrees with what mizerReef
#' does when it projects the model. Predators blocked by refuge
#' (`blocked_pred = TRUE`) only encounter the prey that are not hidden in
#' refuge, `getVulnerable() * n`. mizer's own method does not know about
#' refuge and would compute their diet from the whole prey population,
#' overstating how much they eat of every group that uses refuge. Because
#' [mizer::plotDiet()] calls `getDiet()`, it shows the refuge-aware diet too.
#'
#' The diet of predators that are not blocked by refuge comes from the whole
#' prey population, as in mizer. For every consumer, the part of the
#' encounter that is eaten uses the model's feeding level at the given
#' abundances, algae and detritus biomasses and time.
#'
#' @param object A `mizerReef` params object or a `mizerReefSim` object
#' @param proportion If `TRUE` (default) the function returns the diet as a
#'   proportion of the total consumption rate. If `FALSE` it returns the
#'   consumption rate in grams per year.
#' @param n A matrix of species abundances (species x size). Defaults to the
#'   initial abundances.
#' @param n_pp A vector of the resource abundance by size. Defaults to the
#'   initial resource abundance.
#' @param n_other A list of abundances for other dynamical components, such
#'   as algae and detritus. Defaults to the initial values.
#' @param t For a params object, the time at which the refuge and the feeding
#'   level are calculated. It only matters when refuge degradation is switched
#'   on (see [setDegradation()]). For a `mizerReefSim`, each saved time is
#'   used; choose times with `time_range`.
#' @param time_range The times for which to return the diet, for a
#'   `mizerReefSim`. Defaults to all saved times.
#' @param drop If `TRUE`, dimensions of length 1 are removed from the array
#'   returned for a `mizerReefSim`.
#' @param ... Passed on to [mizer::getDiet()].
#'
#' @return For a params object, an array (predator species x predator size x
#'   prey), as returned by [mizer::getDiet()]. For a `mizerReefSim`, an array
#'   with an additional first dimension for time.
#'
#' @seealso [getVulnerable()]
#' @concept refugeRates
#' @method getDiet mizerReef
#' @export
getDiet.mizerReef <- function(object, proportion = TRUE,
                              n = initialN(object),
                              n_pp = initialNResource(object),
                              n_other = initialNOther(object),
                              t = 0, ...) {
    params <- validParams(object)
    # mizer's method multiplies the encounter by (1 - feeding level), but
    # computes that feeding level from the initial algae and detritus at
    # t = 0. So pass it a copy without satiation, which makes the factor 1,
    # and apply the model's own feeding level below. The model is unchanged.
    feeding_level <- getFeedingLevel(params, n = n, n_pp = n_pp,
                                     n_other = n_other, time_range = t)
    object <- params
    object@intake_max[] <- Inf

    # Name all four arguments. NextMethod() ignores a new value for an
    # argument the caller left out, and naming only some of them shifts
    # positional arguments (see the note in reef-project_methods.R). With all
    # four named, positional copies just end up in the next method's `...`.
    diet <- NextMethod(proportion = FALSE, n = n, n_pp = n_pp,
                       n_other = n_other)

    blocked <- params@species_params$blocked_pred %in% TRUE
    if (any(blocked)) {
        # Blocked predators only encounter the prey outside refuge
        n_vul <- getVulnerable(params, n = n, n_pp = n_pp, n_other = n_other,
                               time_range = t) * n
        # mizer's method zeroes the diet of predator sizes where n is 0, so
        # keep sizes that refuge hides completely just above zero. As prey
        # they then contribute nothing measurable.
        hidden <- n > 0 & n_vul == 0
        n_vul[hidden] <- pmax(n[hidden] * 1e-20, .Machine$double.xmin)
        diet_vul <- NextMethod(proportion = FALSE, n = n_vul, n_pp = n_pp,
                               n_other = n_other)
        fish <- seq_len(nrow(n))
        diet[blocked, , fish] <- diet_vul[blocked, , fish, drop = FALSE]
    }

    diet <- sweep(diet, c(1, 2), 1 - feeding_level, "*")
    if (proportion) {
        total <- rowSums(diet, dims = 2)
        diet <- sweep(diet, c(1, 2), total, "/")
        diet[is.nan(diet)] <- 0
    }
    diet
}

#' @rdname getDiet.mizerReef
#' @method getDiet mizerReefSim
#' @export
getDiet.mizerReefSim <- function(object, proportion = TRUE, time_range,
                                 drop = FALSE, ...) {
    # mizer's MizerSim method does not pass each saved time on, which the
    # refuge and the feeding level need. If mizer passed it, this method
    # could call NextMethod() instead.
    sim <- object
    if (missing(time_range)) {
        time_range <- dimnames(sim@n)$time
    }
    time_elements <- mizer::get_time_elements(sim, time_range)
    diet_time <- plyr::aaply(which(time_elements), 1, function(i) {
        n <- array(sim@n[i, , ], dim = dim(sim@n)[2:3])
        dimnames(n) <- dimnames(sim@n)[2:3]
        n_other <- sim@n_other[i, ]
        names(n_other) <- dimnames(sim@n_other)$component
        getDiet(sim@params, proportion = proportion, n = n,
                n_pp = sim@n_pp[i, ], n_other = n_other,
                t = as.numeric(dimnames(sim@n)$time[[i]]), ...)
    }, .drop = FALSE)
    names(dimnames(diet_time))[[1]] <- "time"
    diet_time[, , , , drop = drop]
}

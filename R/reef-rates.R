#' Get vulnerability level
#'
#' Returns the proportion of fish at size \eqn{w} that are not hidden in
#' predation refuge and thus vulnerable to being encountered by predators.
#'
#' @inherit reefVulnerable
#'
#' @export
#' @family rate functions
getVulnerable <- function(params, n = initialN(params),
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

    f <- get(params@rates_funcs$Vulnerable)
    vulnerable <- f(params, n = n, n_pp = n_pp, n_other = n_other, t = t)
    dimnames(vulnerable) <- dimnames(params@metab)
    vulnerable
}

#' Get encounter rate
#'
#' Returns the rate at which a predator of species \eqn{i} and
#' weight \eqn{w} encounters food (grams/year).
#'
#' @inherit reefEncounter
#' @export
#' @family rate functions
getReefEncounter <- function(params, n = initialN(params),
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
    f <- get(params@rates_funcs$Encounter)
    encounter <- f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                   vulnerable = getVulnerable(params, n = n, n_pp = n_pp,
                                              n_other = n_other,
                                              time_range = t))
    dimnames(encounter) <- dimnames(params@metab)
    encounter
}

# environment(getReefEncounter) <- asNamespace("mizer")
# ## Attempts to fix function replacement issue
# #utils::assignInNamespace("getEncounter", getReefEncounter, ns = "mizer")
# #utils::assignInNamespace("getEncounter", getReefEncounter, pos = "package:mizer")
# #utils::fixInNamespace("getEncounter",  ns = "mizer")
# unlockBinding("getEncounter", as.environment("mizer"))
# assignInNamespace("getEncounter", getReefEncounter, ns= "mizer",
#                   envir=as.environment("mizer"))
# assign("getEncounter", getReefEncounter, as.environment("mizer"))
# lockBinding("getEncounter", as.environment("mizer"))


#' Get predation rate
#'
#' Calculates the potential rate (in units 1/year) at which a prey individual of
#' a given size \eqn{w} is killed by predators from species \eqn{j}. In formulas
#' \deqn{{\tt pred\_rate}_j(w_p) = \int \phi_j(w,w_p) (1-f_j(w))
#'   \gamma_j(w) N_j(w) \, dw.}{pred_rate_j(w_p) = \int\phi_i(w,w_p) (1-f_i(w))
#'   \gamma_i(w) N_i(w) dw.}
#' This potential rate is used in [getPredMort()] to
#' calculate the realised predation mortality rate on the prey individual.
#'
#' @inherit reefPredRate
#'
#' @return A two dimensional array (predator species x prey size),
#'   where the prey size runs over fish community plus resource spectrum.
#' @export
#' @family rate functions
getReefPredRate <- function(params, n = initialN(params),
                            n_pp = initialNResource(params),
                            n_other = initialNOther(params),
                            t = 0, ...) {
    params <- validParams(params)
    f <- get(params@rates_funcs$PredRate)
    pred_rate <- f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
                   feeding_level = getFeedingLevel(params, n = n, n_pp = n_pp,
                                                   n_other = n_other,
                                                   time_range = t),
                   vulnerable = getVulnerable(params, n = n, n_pp = n_pp,
                                              n_other = n_other,
                                              time_range = t))
    dimnames(pred_rate) <- list(sp = dimnames(params@initial_n)$sp,
                                w_prey = as.character(signif(params@w_full, 3)))
    pred_rate
}

# environment(getReefPredRate) <- asNamespace("mizer")
# #utils::assignInNamespace("getPredRate", getReefPredRate, ns = "mizer")
# #utils::assignInNamespace("getPredRate", getReefPredRate, pos = "package:mizer")
# #utils::fixInNamespace("getPredRate",  ns = "mizer")
# unlockBinding("getPredRate", as.environment("mizer"))
# assignInNamespace("getPredRate", getReefPredRate, ns= "mizer",
#                   envir=as.environment("mizer"))
# assign("getPredRate", getReefPredRate, as.environment("mizer"))
# lockBinding("getPredRate", as.environment("mizer"))

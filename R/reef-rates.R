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
                is.number(t),
                identical(dim(n), dim(params@initial_n)),
                identical(length(n_pp), length(params@initial_n_pp)),
                identical(length(n_other), length(params@initial_n_other))
    )

    vulnerable <- reefVulnerable(params, n = n, n_pp = n_pp, n_other = n_other, t = t)
    dimnames(vulnerable) <- dimnames(params@metab)
    vulnerable
}



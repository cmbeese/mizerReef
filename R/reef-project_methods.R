#' Get all rates needed to project a mizerReef model
#'
#' Calls other rate functions in sequence and collects the results in a list.
#'
#' By default this function returns a list with the following components:
#'
#'  \itemize{
#'      \item predation vulnerability from [reefVulnerable()]
#'      \item encounter from [reefEncounter()]
#'      \item feeding level from [reefFeedingLevel()]
#'      \item e from [mizerEReproAndGrowth()]
#'      \item e_repro from [mizerERepro()]
#'      \item e_growth from [mizerEGrowth()]
#'      \item pred_rate from [reefPredRate()]
#'      \item pred_mort from [mizerPredMort()]
#'      \item sen_mort from [reefSenMort()]
#'      \item f_mort from [mizerFMort()]
#'      \item mort from [reefMort()]
#'      \item rdi from [mizerRDI()]
#'      \item rdd from [BevertonHoltRDD()]
#'      \item resource_mort from [mizerResourceMort()]
#'   }
#'
#' However you can replace any of these rate functions by your own rate
#' function if you wish, see [setRateFunction()] for details.
#'
#' @param params A \linkS4class{MizerParams} object
#' @param n A matrix of species abundances (species x size).
#' @param n_pp A vector of the resource abundance by size
#' @param n_other A list of abundances for other dynamical components of the
#'   ecosystem
#' @param t The time for which to do the calculation (Not used by standard
#'   mizer rate functions but useful for extensions with time-dependent
#'   parameters.)
#' @param effort The effort for each fishing gear
#' @param rates_fns Named list of the functions to call to calculate the rates.
#'   Note that this list holds the functions themselves, not their names.
#' @param ... Unused
#' @return List of rates.
#' @export
#' @concept refuge
#' @family mizer rate functions
reefRates <- function(params, n, n_pp, n_other,
                      t = 0, effort, rates_fns, ...) {
    r <- list()

    ## Vulnerability ----
    # Calculate vulnerability of fish based on complexity
    r$vulnerable <- reefVulnerable(
        params, n = n, n_pp = n_pp, n_other = n_other, t = t, ...)
    
    # Implement degradation
    r$degrade <- reefDegrade(
        params, n = n, n_pp = n_pp, n_other = n_other, 
        vulnerable = r$vulnerable, t = t, ...)
    
    ## Growth ----
    # Calculate rate E_{e,i}(w) of encountered food
    r$encounter <- rates_fns$Encounter(
        params, n = n, n_pp = n_pp, n_other = n_other,
        vulnerable = r$vulnerable, t = t, ...)
    # Calculate feeding level f_i(w)
    r$feeding_level <- rates_fns$FeedingLevel(
        params, n = n, n_pp = n_pp, n_other = n_other,
        encounter = r$encounter, t = t, ...)
    # Calculate the energy available for reproduction and growth
    r$e <- rates_fns$EReproAndGrowth(
        params, n = n, n_pp = n_pp, n_other = n_other,
        encounter = r$encounter, feeding_level = r$feeding_level, t = t, ...)
    # Calculate the energy for reproduction
    r$e_repro <- rates_fns$ERepro(
        params, n = n, n_pp = n_pp, n_other = n_other,
        e = r$e, t = t, ...)
    # Calculate the growth rate g_i(w)
    r$e_growth <- rates_fns$EGrowth(
        params, n = n, n_pp = n_pp, n_other = n_other,
        e_repro = r$e_repro, e = r$e, t = t, ...)

    ## Mortality ----
    # Calculate the predation rate
    r$pred_rate <- rates_fns$PredRate(
        params, n = n, n_pp = n_pp, n_other = n_other,
        feeding_level = r$feeding_level, vulnerable = r$vulnerable,
        t = t, ...)
    # Calculate predation mortality on fish \mu_{p,i}(w)
    r$pred_mort <- rates_fns$PredMort(
        params, n = n, n_pp = n_pp, n_other = n_other,
        pred_rate = r$pred_rate, t = t, ...)
    # Calculate fishing mortality
    r$f_mort <- rates_fns$FMort(
        params, n = n, n_pp = n_pp, n_other = n_other,
        effort = effort, t = t,
        e_growth = r$e_growth, pred_mort = r$pred_mort, ...)
    # Calculate total mortality \mu_i(w)
    r$mort <- rates_fns$Mort(
        params, n = n, n_pp = n_pp, n_other = n_other,
        f_mort = r$f_mort, pred_mort = r$pred_mort, t = t, ...)

    ## Reproduction ----
    # R_di
    r$rdi <- rates_fns$RDI(
        params, n = n, n_pp = n_pp, n_other = n_other,
        e_growth = r$e_growth,
        mort = r$mort,
        e_repro = r$e_repro, t = t, ...)
    # R_dd
    r$rdd <- rates_fns$RDD(
        rdi = r$rdi, species_params = params@species_params, ...)

    ## Resource ----
    # Calculate mortality on the resource spectrum
    r$resource_mort <- rates_fns$ResourceMort(
        params, n = n, n_pp = n_pp, n_other = n_other,
        pred_rate = r$pred_rate, t = t, ...)

    return(r)
}

#' Define proportion of fish hidden in predation refuge to simulate 
#' benthic habitat complexity
#'
#' This function calculates the proportion of fish that are not hidden in
#' predation refuge and thus vulnerable to being encountered by predators.
#' 
#' @inheritSection setRefuge Setting the refuge profile'
#' @inheritParams reefRates
#' @param ... Unused
#'
#' @return Array (species x size) with the proportion of individuals that are
#'          not protected from predation by refuge
#'
#' @export
#' @concept refuge
#' @family mizer rate functions
reefVulnerable <- function(params, n, n_pp, n_other, t = 0, ...) {

    # Extract relevant data from params
    refuge_params <- params@other_params[['refuge_params']]
    method_params <- params@other_params[['method_params']]

    # Pull values from params
    w <- params@w
    no_w <- length(params@w)
    no_sp <- dim(params@interaction)[1]

    # Initialize storage for the array of refuge proportions
    refuge <- matrix(0, nrow = no_sp, ncol = no_w)
    rownames(refuge) <- rownames(params@initial_n)
    colnames(refuge) <- colnames(params@initial_n)

    # Set parameters used with all methods
    w_settle    <- refuge_params$w_settle
    # w_settle    <- params@species_params$w_settle
    max_protect  <- refuge_params$max_protect
    tau          <- refuge_params$tau

    # Store which functional groups use refuge
    refuge_user <- params@species_params$refuge_user

    # Static methods
    # static = c("sigmoidal", "binned")
    # if (is.element(refuge_params$method, static)){
    #
    #     vulnerable <- params@other_params$initial_vulnerable

    # Sigmoidal method ------------------------------------------------------------
    if (refuge_params$method == "sigmoidal"){

        # Pull slop and proportion of fish to be protected from method_params
        prop_protect <- method_params$prop_protect
        slope <- method_params$slope

        # Convert length to weight
        W_refuge <- params@species_params[["a"]] *
            method_params$L_refuge ^ params@species_params[["b"]]

        # Find indices of fish in size range to protect
        # Set threshold weight - no organisms smaller than min_ref_length
        # can utilize refuge to escape predators
        # Loop through species
        for (i in 1:no_sp){
            denom <- 1 + exp(slope*(w - W_refuge[i]))
            refuge[i, ] <- ifelse(w > w_settle, prop_protect/denom, 0)
        }

    # Binned method ------------------------------------------------------------
    } else if (refuge_params$method == "binned") {

        # Loop through each refuge bin
        for (k in 1:nrow(method_params)) {

            # Calculate start and end of weight bins for each functional group
            # based on unique as and bs
            start_w <- params@species_params[["a"]] *
                method_params$start_L[[k]] ^ params@species_params[["b"]]

            end_w  <- params@species_params[["a"]] *
                method_params$end_L[[k]] ^ params@species_params[["b"]]

            # Set threshold weight - no organisms smaller than w_settle
            # can utilize refuge to escape predators
            # idx.sm <- which(start_w < w_settle)
            # start_w[idx.sm] <- w_settle[idx.sm]
            start_w[start_w < w_settle] <- w_settle

            # Find indices of fish in size range to protect
            # TRUE for fish larger than start weight of bin
            start_n <- t(sapply(start_w, function(x) params@w >= x))
            # TRUE for fish smaller than end weight of bin
            end_n <- t(sapply(end_w, function(x) params@w <= x))
            # TRUE for fish that meet both conditions
            bin_fish <- start_n & end_n

            # Find indices of functional group x weight that are in size bin k
            idx.bin = which(bin_fish == TRUE)

            # Set vulnerability for fish in size bin to provided value
            refuge[idx.bin] = method_params$prop_protect[k]
        }

    # Data method --------------------------------------------------------------
    } else if (refuge_params$method == "data") {

        # Initialize empty list to hold number of competitors for each bin
        competitor_density = numeric(nrow(method_params))

        # Loop through each refuge bin
        for (k in 1:nrow(method_params)) {

            # Calculate start and end of weight bins for each functional
            # group based on unique as and bs
            start_w <- params@species_params[["a"]] *
                method_params$start_L[[k]] ^ params@species_params[["b"]]

            end_w  <- params@species_params[["a"]] *
                method_params$end_L[[k]] ^ params@species_params[["b"]]

            # Set threshold weight - no organisms smaller than
            # minimum size that can utilize refuge
            # idx.sm <- which(start_w < w_settle)
            # start_w[idx.sm] <- w_settle[idx.sm]
            start_w[start_w < w_settle] <- w_settle

            # Calculate competitor density - number of fish that use refuge in
            # each size bin

            # Create matrix of boolean values indicating whether fish is
            # in size range of bin k
            # TRUE for fish larger than start weight
            start_n <- t(sapply(start_w, function(x) params@w >= x))
            # TRUE for fish smaller than end weight
            end_n <- t(sapply(end_w, function(x) params@w <= x))
            # TRUE for fish that meet both conditions
            bin_fish <- start_n & end_n

            # Number of competitors from each functional group in bin
            competitors <- (n * bin_fish) %*% params@dw

            # Eliminate functional groups that don't use refuge and sum
            competitor_density[k] <- sum(refuge_user * competitors)

            # Find indices of fish within size bin k
            ind = which(bin_fish == TRUE)

            # Set vulnerability for fish in size bin based on the number of
            # available refuges and the number of competitors
            refuge[ind] <- ifelse(competitor_density[k] == 0, max_protect,
                                  tau * method_params$refuge_density[k] /
                                      competitor_density[k])
        }
    }

    # Save results
    # Make sure none of the values are higher than maximum protection allowed
    refuge[refuge > max_protect] = max_protect

    # Account for species that don't utilize refuge
    vulnerable = 1 - (refuge_user*refuge)

    return(vulnerable)
}

# #' Degrade coral reef habitat structure by decreasing the availability of 
# #' refuge at certain sizes
# #'
# #' TO DO: ADD CODE AND DESCRIPTION HERE OF DEGRADATION
# #'
# #' @inheritParams reefRates
# #' @param vulnerable A two dimensional array (prey species x prey size) with
# #'      the proportion of prey vulnerable to being encountered.
# #' @param ... Unused
# #'
# #' @return Array (species x size) with the proportion of individuals that are
# #'          not protected from predation by refuge
# #'
# #' @export
# #' @family mizer rate functions
# # Account for species that don't utilize refuge
# reefDegrade <- function(params, n, n_pp, n_other, t, 
#                         vulnerable = reefVulnerable(params, 
#                                                     n, n_pp, n_other, t),...)
#     {
#     
#     if (degrade == FALSE){
#         return(vulnerable)
#     } else {
#         
#         decrease_amount <- (percentage / 100) * value
#         result <- value - decrease_amount
#         
#     }   
#     
# }


#' Get encounter rate needed to project a mizerReef model
#'
#' Calculates the rate \eqn{E_i(w)} at which a predator from group \eqn{i} and
#' weight \eqn{w} encounters food (grams/year). You would not usually call this
#' function directly but instead use [getEncounter()], which then calls this
#' function.
#'
#' @section Predation encounter:
#'
#'      The encounter rate \eqn{E_i(w)} at which a predator of species \eqn{i}
#'      and weight \eqn{w} encounters food has contributions from the encounter
#'      of fish prey and of resources. This is determined by summing over all
#'      prey species and the resource spectrum and then integrating over all
#'      prey sizes \eqn{w_p}, weighted by predation kernel \eqn{\phi(w,w_p)}:
#'
#'           \deqn{
#'          E_i(w) = \gamma_i(w) \int
#'          \left( \theta_{ip} N_R(w_p) + \sum_{j} \theta_{ij}
#'          V_{ji}(w_p) N_j(w_p) \right)
#'          \phi_i(w,w_p) w_p \, dw_p.
#'          }{\gamma_i(w) \int
#'          ( \theta_{ip} N_R(w_p) + \sum_{j} \theta_{ij} V_{ji}(w_p) N_j(w_p) )
#'          \phi_i(w,w_p) w_p dw_p.}
#'
#'      Here \eqn{N_j(w)} is the abundance density of species \eqn{j} and
#'      \eqn{N_R(w)} is the abundance density of resource. The overall
#'      prefactor \eqn{\gamma_i(w)} determines the predation power of the
#'      predator. It could be interpreted as a search volume and is set with
#'      the [setSearchVolume()] function.
#'
#'      The predation kernel \eqn{\phi(w,w_p)}is set with the [setPredKernel()]
#'      function.
#'
#'      The vulnerability to predation, \eqn{V_{ji}(w)} accounts for protective
#'      behavior of the prey. The parameters that control this are set with the
#'      [setRefuge()] function.
#'
#'      The species interaction matrix \eqn{\theta_{ij}} is set with
#'      [setInteraction()] and the resource interaction vector \eqn{\theta_{ip}}
#'      is taken from the `interaction_resource`column in
#'      `params@species_params`.
#'
#' @section Details:
#' The encounter rate is multiplied by \eqn{1-f_0} to obtain the consumption
#' rate, where \eqn{f_0} is the feeding level calculated with
#' [getFeedingLevel()]. This is used by the [project()] function for performing
#' simulations.
#'
#' The function returns values also for sizes outside the size-range of the
#' species. These values should not be used, as they are meaningless.
#'
#' If your model contains additional components that you added with
#' [setComponent()] and for which you specified an `encounter_fun` function then
#' the encounters of these components will be included in the returned value.
#'
#' @inheritParams reefRates
#' @param vulnerable A two dimensional array (prey species x prey size) with
#'      the proportion of prey vulnerable to being encountered.
#' @param ... Unused
#'
#' @return A named two dimensional array (predator species x predator size) with
#'   the encounter rates.
#' @export
#' @concept refuge
#' @family mizer rate functions
reefEncounter <- function(params, n, n_pp, n_other, t, 
                          vulnerable = reefVulnerable(params, 
                                                      n, n_pp, n_other, t),...) {

    # Pull values from params
    no_sp <- dim(params@interaction)[1]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)

    # idx_sp are the index values of params@w_full such that
    # params@w_full[idx_sp] = params@w
    idx_sp <- (no_w_full - no_w + 1):no_w_full

    # Initialize encounter matrix
    enc <- matrix(0, no_sp, no_w)

    # Find indices of predator species impacted by refuge
    bad_pred  <- which(params@species_params$bad_pred == TRUE)
    good_pred <- which(params@species_params$bad_pred == FALSE)

    # Calculate n_vulnerable, number at each size vulnerable to being
    # encountered
    n_vul <- vulnerable * n

    # If the the user has set a custom pred_kernel we can not use fft.
    # In this case we use the code from mizer version 0.3
    if (!is.null(comment(params@pred_kernel))) {
        # n_eff_prey is the total prey abundance by size exposed to each
        # predator (prey not broken into species - here we are just working out
        # how much a predator eats - not which species are being eaten - that is
        # in the mortality calculation
        # \sum_j \theta_{ij} N_j(w_p) w_p dw_p

        # First deal with encounter rate for predators unaffected by refuge
        n_eff_prey <- sweep(params@interaction %*% n, 2,
                            params@w * params@dw, "*", check.margin = FALSE)
        # pred_kernel is predator species x predator size x prey size
        # So multiply 3rd dimension of pred_kernel by the prey biomass density
        # Then sum over 3rd dimension to get consumption rate of each predator
        # by predator size
        # This line is a bottle neck
        phi_prey_species <- rowSums(sweep(
            params@pred_kernel[, , idx_sp, drop = FALSE],
            c(1, 3), n_eff_prey, "*", check.margin = FALSE), dims = 2)
        # Eating the background
        # This line is a bottle neck
        phi_prey_background <- params@species_params$interaction_resource *
            rowSums(sweep(
                params@pred_kernel, 3, params@dw_full * params@w_full * n_pp,
                "*", check.margin = FALSE), dims = 2)

        good_encounter <- params@search_vol * (phi_prey_species +
                                                   phi_prey_background)
        enc[good_pred,] <- good_encounter[good_pred,]

        # Now deal with predators who are affected by refuge
        v_n_eff_prey <- sweep(params@interaction %*% n_vul, 2,
                              params@w * params@dw, "*", check.margin = FALSE)

        v_phi_prey_species <- rowSums(sweep(
            params@pred_kernel[, , idx_sp, drop = FALSE],
            c(1, 3), v_n_eff_prey, "*", check.margin = FALSE), dims = 2)

        bad_encounter <- params@search_vol * (v_phi_prey_species +
                                                  phi_prey_background)

        enc[bad_pred,] <- bad_encounter[bad_pred,]

        encounter <- enc

    } else {

        # First deal with predators that are not affected by refuge
        prey <- outer(params@species_params$interaction_resource, n_pp)
        prey[, idx_sp] <- prey[, idx_sp] + params@interaction %*% n
        # The vector prey equals everything inside integral (3.4)
        # except the feeding kernel phi_i(w_p/w).
        prey <- sweep(prey, 2, params@w_full * params@dw_full, "*")
        # Eq (3.4) is then a convolution integral in terms of prey[w_p]
        # and phi[w_p/w]. We approximate the integral by the trapezoidal
        # method. Using the convolution theorem we can evaluate the resulting
        # sum via fast fourier transform.
        # mvfft() does a Fourier transform of each column of its argument, but
        # we need the Fourier transforms of each row, so we need to apply
        # mvfft() to the transposed matrices and then transpose again at
        # the end.
        avail_energy <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_e) *
                                             mvfft(base::t(prey)),
                                    inverse = TRUE))) / length(params@w_full)
        # Only keep the bit for fish sizes
        avail_energy <- avail_energy[, idx_sp, drop = FALSE]
        # Due to numerical errors we might get negative or very small entries
        # that should be 0
        avail_energy[avail_energy < 1e-18] <- 0

        good_enc <- params@search_vol * avail_energy
        enc[good_pred,] <- good_enc[good_pred,]

        # Now deal with predators who are affected by refuge
        v_prey <- outer(params@species_params$interaction_resource, n_pp)
        v_prey[, idx_sp] <- v_prey[, idx_sp] + params@interaction %*% n_vul
        # The vector prey equals everything inside integral (3.4)
        # except the feeding kernel phi_i(w_p/w).
        v_prey <- sweep(v_prey, 2, params@w_full * params@dw_full, "*")

        v_avail_energy <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_e) *
                                               mvfft(base::t(v_prey)),
                                           inverse = TRUE)))/length(params@w_full)

        # Only keep the bit for fish sizes
        v_avail_energy <- v_avail_energy[, idx_sp, drop = FALSE]
        # Due to numerical errors we might get negative or very small entries
        # that should be 0
        v_avail_energy[v_avail_energy < 1e-18] <- 0

        bad_enc <- params@search_vol * v_avail_energy
        enc[bad_pred,] <- bad_enc[bad_pred,]

        encounter <- enc
    }

    # Add contributions from other components
    for (i in seq_along(params@other_encounter)) {
        encounter <- encounter +
            do.call(params@other_encounter[[i]],
                    list(params = params,
                         n = n, n_pp = n_pp, n_other = n_other,
                         component = names(params@other_encounter)[[i]], ...))
    }

    return(encounter)
}

#' Get total predation mortality rate needed to project mizer reef model
#'
#' Calculates the total predation mortality rate \eqn{\mu_{p,i}(w_p)} (in units
#' of 1/year) on each prey species by prey size:
#' 
#' \deqn{\mu_{p.i}(w_p) = \sum_j {\tt pred\_rate}_j(w_p)\, V_{ji}(w_p)\, \theta_{ji}.}{
#'   \mu_{p.i}(w_p) = \sum_j pred_rate_j(w_p) V_{ji}(w_p) \theta_{ji}.}
#'   
#' You would not usually call this
#' function directly but instead use [getPredMort()], which then calls this
#' function unless an alternative function has been registered, see below.
#' 
#' @inheritParams reefRates
#' @param pred_rate A two dimensional array (predator species x predator size)
#'   with the feeding level.
#'
#' @return A two dimensional array (prey species x prey size) with the predation
#'   mortality
#' @family mizer rate functions
#' @concept refuge
#' @export
reefPredMort <- function(params, n, n_pp, n_other, t, pred_rate,
                         vulnerable = reefVulnerable(params, n, n_pp,
                                                     n_other, t), ...) {
    
    # Find indices of fish that have grown out of the resource spectrum
    idx_sp <- (length(params@w_full) - 
                   length(params@w) + 1):length(params@w_full)
    no_sp <- nrow(params@species_params)
    
    # Find indices of predator species whose foraging is hindered by refuge
    bad_pred  <- which(params@species_params$bad_pred == TRUE)
    good_pred <- which(params@species_params$bad_pred == FALSE)
    
    # Create list of vulnerabilities for each predator
    vul <- vector("list", no_sp)
    vul[bad_pred] <- list(vulnerable)
    vul[good_pred] <- list(1)
    
    # Loop through predator species to calculate predation mortality on
    # each prey species & size by all predators
    for (i in 1:no_sp){
        # Predation rate of predator species i on all prey (by size)
        pr_i <- pred_rate[i, idx_sp, drop = TRUE]
        # Vulnerability rate of all prey (species by size) to predator i
        vul <- vul[[i]]
        # vul*pr_i predation mortality on prey (species by size) by predator i 
        # pm + add to predation mortality from other predators
        pm <- pm + vul*pr_i
    }
    
    return(base::t(params@interaction) %*% pm)
}


#' Expanding external mortality rate to include senescence
#'
#' @section Senescence mortality:
#'
#'      Senescence mortality \eqn{\mu_{sen.i}(w)} is used to represent
#'      mortality caused by external sources such as illness or age. This is
#'      addition to external mortality, \eqn{\mu_{ext.i}(w)}, which represents
#'      all mortality that is not due to fishing or predation by predators
#'      included in the model. The rate of senescence mortality is given by:
#'
#'      \deqn{\mu_{sen.i}(w) = k_{sen}\left(\frac{w}{w_{max.i}}\right)^{p_{sen}}}
#'           {\mu_{sen.i}(w) = k_{sen} (w/w_{max.i})^{p_{sen}}
#'
#'      where \eqn{k_{sen}} and \eqn{p_sen} are constants defining the shape
#'      of the senescence curve and \eqn{w_max.i} is the maximum size
#'      of species \eqn{i} in grams.
#'
#'      Users can change all constants with the `setSenMortParams()` function.
#'
#' @param params A MizerParams object
#' @param n A matrix of species abundances (species x size).
#' @param n_pp A vector of the resource abundance by size
#' @param n_other A list of abundances for other dynamical components of the
#'   ecosystem
#' @param t The time for which to do the calculation (Not used by standard
#'   mizer rate functions but useful for extensions with time-dependent
#'   parameters.)
#'
#' @return A named two dimensional array (species x size) with the senescence
#'   mortality rates.
#' @concept External mortality
#' @export
reefSenMort <- function(params, ...) {

    # Pull values from params for use later
    no_sp <- dim(params@interaction)[1]
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)

    # Get user set senescence mortality parameters
    mort_params <- params@other_params[['ext_mort_params']]
    sen_prop    <- mort_params$sen_prop
    sen_curve   <- mort_params$sen_curve

    # Initialize storage for each mortality rates
    sen_mort <- matrix(0, nrow = no_sp, ncol = no_w)

    # # Convert length to weight to get senesence weight for each species
    # sen_weight <- params@species_params[["a"]] *
    #     sen_length ^ params@species_params[["b"]]

    # Or could use max size?
    sen_weight <- params@species_params[['w_max']]
    log_sw <- log10(sen_weight)
    log_w  <- log10(params@w)

    # Loop through species
    for (i in 1:length(sen_weight)){
        sen <- sen_prop * ( log_w / log_sw[i] )
        sen[sen < 0] <- 0 
        sen_mort[i,] <- sen ^ sen_curve
    }
    return(sen_mort)
}

#' Total mortality rate in the reef ecosystem model
#'
#' This function replaces the usual [mizerMort()] function and returns the
#' sum of the usual mortality and size-based external/ senescence mortality
#'
#' @param params A MizerParams object
#' @param n A matrix of species abundances (species x size).
#' @param n_pp A vector of the resource abundance by size
#' @param n_other A list of abundances for other dynamical components of the
#'   ecosystem
#' @param t The time for which to do the calculation (Not used by standard
#'   mizer rate functions but useful for extensions with time-dependent
#'   parameters.)
#' @param f_mort A two dimensional array (species x size) with the fishing
#'   mortality
#' @param pred_mort A two dimensional array (species x size) with the predation
#'   mortality
#' @param ... Unused
#'
#' @return A named two dimensional array (species x size) with the total
#'   mortality rates.
#' @family mizer rate functions
#' @concept External mortality
#' @export
reefMort <- function(params, n, n_pp, n_other, t, f_mort, pred_mort, ...) {
    mizerMort(params, n, n_pp, n_other, t, f_mort, pred_mort, ...) +
        reefSenMort(params, ...)
}




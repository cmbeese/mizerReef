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
#'      \item pred_rate from [mizerPredRate()]
#'      \item pred_mort from [reefPredMort()]
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
#' 
#' @param n A matrix of species abundances (species x size).
#' 
#' @param n_pp A vector of the resource abundance by size
#' 
#' @param n_other   A list of abundances for other dynamical components of the
#'                  ecosystem
#'                  
#' @param t The time for which to do the calculation (Not used by standard
#'          mizer rate functions but useful for extensions with time-dependent
#'          parameters.)
#'          
#' @param effort The effort for each fishing gear
#' 
#' @param rates_fns Named list of the functions to call to calculate the rates.
#'                  Note that this list holds the functions themselves, not 
#'                  their names.
#'                  
#' @param ... Unused
#' @return List of rates.
#' @export
#' @concept refugeRates
#' @family mizer rate functions
reefRates <- function(params, n, n_pp, n_other,
                      t = 0, effort, rates_fns, ...) {
    r <- list()
    
    dt <- params@other_params$dt
    ## Degradation -----
    # Implement degradation in projections
    r$degrade <- reefDegrade(
        params, n = n, n_pp = n_pp, n_other = n_other, t = t,
        old_rd = reefDegrade(params, n, n_pp, n_other,
                             t = t, dt = dt), ...)
    
    ## Vulnerability ----
    # Calculate vulnerability of fish based on complexity
    r$vulnerable <- reefVulnerable(
        params, n = n, n_pp = n_pp, n_other = n_other, t = t, 
        new_rd = r$degrade, ...)
    
    ## Growth ----
    # Calculate rate E_{e,i}(w) of encountered food
    r$encounter <- reefEncounter(
        # rates_fns$Encounter(
        params, n = n, n_pp = n_pp, n_other = n_other,
        vulnerable = r$vulnerable, t = t, ...)
    # Calculate feeding level f_i(w)
    r$feeding_level <- reefFeedingLevel(
        # rates_fns$FeedingLevel(
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
    r$pred_mort <- reefPredMort(
        # rates_fns$PredMort(
        params, n = n, n_pp = n_pp, n_other = n_other,
        pred_rate = r$pred_rate, t = t, ...)
    # Calculate fishing mortality
    r$f_mort <- rates_fns$FMort(
        params, n = n, n_pp = n_pp, n_other = n_other,
        effort = effort, t = t,
        e_growth = r$e_growth, pred_mort = r$pred_mort, ...)
    # Calculate total mortality \mu_i(w)
    r$mort <- reefMort(
        # rates_fns$Mort(
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

#' Scales the refuge density by a given value at set times
#' 
#' Allows for the degradation of coral reef habitat structure following an
#' acute disturbance by decreasing the availability of refuge over time. 
#'
#' @inheritParams reefRates
#' @param ... Unused
#'
#' @return A new methods parameters data frame scaled by bleaching
#' @concept degradation
#' @export
#' @family mizer rate functions
# Check with Gustav about hard coding dt, see if its stored anywhere or can
# be use in rate functions
reefDegrade <- function(params, n, n_pp, n_other, t,...) {
    
    method_params <- params@other_params[['method_params']]
    refuge_params <- params@other_params[['refuge_params']]
    degrade <- params@other_params$degrade
    
    # Save original refuge density
    rd <- method_params$refuge_density
    
        
    if (degrade == TRUE){
        # If not competitive method, error
        if (refuge_params$method != "competitive"){
            stop("Degradation is only available for the competitive
                 method.")
        }
    
        # Pull time for bleaching & scaling parameters
        t_bleach    <- params@other_params$t_bleach
        trajectory  <- params@other_params$trajectory
        deg_scale   <- params@other_params[['deg_scale']]
        
        # If before bleaching, return old method parameters
        time <- t
        if(time < t_bleach) {return(rd)}
        
        old_rd <- reefDegrade(params, n, n_pp, n_other, t = t-1)
        
        # Bleaching time - at bleach year
        if(time == t_bleach) {
            # Set bleach time scale factors
            sev_bleach  <- c(0.1, 0.1, rep(0.9, 8))
            mild_bleach <- c(0.3, 0.7, rep(0.9, 8))
            # Set trajectories that are severe and mild
            sev_traj  <- c("rubble")
            mild_traj <- c("algae", "recovery")
            
            # Scale refuge
            if(trajectory %in% sev_traj) {
                new_rd <-  sev_bleach * old_rd
            }
            if(trajectory %in% mild_traj){
                new_rd <-  mild_bleach * old_rd
            }
            
            # # Algal dynamics
            if (trajectory == "algae"){
                a_growth <- params@other_params$algae$growth
                a_capacity <- params@other_params$algae$capacity
                new_a_growth <- 1.11*a_growth
                new_a_capacity <- 2*a_capacity
                params@other_params$algae$growth <- 1.11*new_a_growth
                params@other_params$algae$capacity <- new_a_capacity
            }
            if(trajectory == "rubble"){
                a_capacity <- params@other_params$algae$carry
                new_a_capacity <- 1.5*a_capacity
                params@other_params$algae$capacity <- new_a_capacity
            }

            # Save new method_params in parameter spot so that scaling
            # occurs on new values next time step
            params@time_modified <- lubridate::now()

            return(new_rd)
        }
        
        if(time > t_bleach){
            # Calculate number of years post bleaching
            years_post <- time - t_bleach
            
            if(years_post %in% 1:dim(deg_scale)[2]) {
                # Get refuge density scaling from deg_scale
                # if in bleaching time range, otherwise scale by 1
                scale_bin <- deg_scale[,years_post]
                # multiply bins by scaling values
                new_rd <- scale_bin * old_rd
                # Increase algae growth and capacity if algae method
                # and in first three years following bleaching
                if (trajectory == "algae"){
                    if (years_post %in% (1:3)) {
                        a_growth <- params@other_params$algae$growth
                        new_a_growth <- 1.11*a_growth
                        params@other_params$algae$growth <- new_a_growth
                    }
                }
                return(new_rd)
            }
            return(old_rd)
        }
    } else {
        return(rd)
    }
}

#' Find the proportion of fish vulnerable to being encountered by predators 
#' at each time step
#'
#' This function calculates the proportion of fish that are not hidden in
#' predation refuge and thus vulnerable to being encountered by predators.
#' 
#' @inheritSection setRefuge Setting the refuge profile
#' @inheritParams reefRates
#' @param ... Unused
#'
#' @return Array (species x size) with the proportion of individuals that are
#'          not protected from predation by refuge
#'
#' @export
#' @concept refugeRates
#' @family mizer rate functions
#' 
reefVulnerable <- function(params, n, n_pp, n_other, t,
    new_rd = reefDegrade(params,  n, n_pp, n_other, t,...)) {
    
    # Extract relevant data from params
    refuge_params <- params@other_params[['refuge_params']]
    method_params <- params@other_params[['method_params']]
    
    method_params$refuge_density <- new_rd
    
    # Set parameters used with all methods
    w_settle    <- refuge_params$w_settle
    max_protect <- refuge_params$max_protect
    tau         <- refuge_params$tau
    
    # Pull no of species and size bins
    no_w <- length(params@w)
    no_sp <- dim(params@interaction)[1]
    
    # Store which functional groups use refuge
    refuge_user <- params@species_params$refuge_user
    
    # Static methods -----------------------------------------------------------
    static = c("sigmoidal", "binned", "noncomplex")
    
    if (is.element(refuge_params$method, static)){
        
        refuge <- params@other_params$refuge
        vulnerable <- 1 - refuge
        
    # Competitive method -------------------------------------------------------
    } else if (refuge_params$method == "competitive") {
        
        # Initialize empty list to hold number of competitors for each bin
        competitor_density = numeric(length(method_params$refuge_density))
        
        # Initialize storage for the array of refuge proportions
        refuge <- matrix(0, nrow = no_sp, ncol = no_w)
        rownames(refuge) <- rownames(params@initial_n)
        colnames(refuge) <- colnames(params@initial_n)
        
        # Loop through each refuge bin
        for (k in 1:nrow(method_params)) {
            # Get indices of fish in size bin k
            bin.id <- params@other_params$bin.id[[k]]
            
            # Create logical vector and use to get abundances in size bin
            bin_fish <- 1:no_w %in% bin.id
            bin_fish <- sweep(n, 2, bin_fish, "*") 
            
            # Calculate number of competitors from each species group in bin k
            competitors <- bin_fish %*% params@dw
            
            # Remove species that don't use refuge
            sp <- params@species_params$species
            sp <- sp[params@species_params$refuge_user ==TRUE]
            competitors <- competitors[sp,]
            
            # sum competitors from all species groups for refuge bin k
            competitor_density[k] <- sum(competitors)
            
            # Set vulnerability for fish in size bin based on the number of
            # available refuges and the number of competitors
            refuge[,bin.id] <- ifelse(competitor_density[k] == 0, 
                max_protect,
                tau * method_params$refuge_density[k]/competitor_density[k])
        }
        # Make sure none of the values are higher than max_protect
        refuge[refuge > max_protect] <- max_protect
        # Account for vulnerability of species that don't utilize refuge
        vulnerable <- 1 - (refuge_user*refuge)
    }
    
    return(vulnerable)
}

#' Get encounter rate needed to project a mizerReef model
#'
#' Calculates the rate \eqn{E_i(w)} at which a predator from group \eqn{i} and
#' weight \eqn{w} encounters food (grams/year). You would not usually call this
#' function directly but instead use [getEncounter()], which then calls this
#' function.
#'
#' @section Predation encounter:
#'
#'  The encounter rate \eqn{E_i(w)} at which a predator of species \eqn{i}
#'  and weight \eqn{w} encounters food has contributions from the encounter
#'  of fish prey and of resources. This is determined by summing over all
#'  prey species and the resource spectrum and then integrating over all
#'  prey sizes \eqn{w_p}, weighted by predation kernel \eqn{\phi(w,w_p)}:
#'  
#'   \deqn{
#'  E_i(w) = \gamma_i(w) \int
#'  \left( \theta_{ip} N_R(w_p) + \sum_{j} \theta_{ij}
#'  V_{ji}(w_p) N_j(w_p) \right)
#'  \phi_i(w,w_p) w_p \, dw_p.
#'  }{\gamma_i(w) \int
#'  ( \theta_{ip} N_R(w_p) + \sum_{j} \theta_{ij} V_{ji}(w_p) N_j(w_p) )
#'  \phi_i(w,w_p) w_p dw_p.}
#
#'  Here \eqn{N_j(w)} is the abundance density of species \eqn{j} and
#'  \eqn{N_R(w)} is the abundance density of resource. The overall
#'  prefactor \eqn{\gamma_i(w)} determines the predation power of the
#'  predator. It could be interpreted as a search volume and is set with
#'  the [setSearchVolume()] function.
#'  
#'  The predation kernel \eqn{\phi(w,w_p)}is set with the [setPredKernel()]
#'  function.
#'  
#'  The vulnerability to predation, \eqn{V_{ji}(w)} accounts for protective
#'  behavior of the prey. The parameters that control this are set with the
#'  [setRefuge()] function.
#'  
#'  The species interaction matrix \eqn{\theta_{ij}} is set with
#'  [setInteraction()] and the resource interaction vector \eqn{\theta_{ip}}
#'  is taken from the `interaction_resource`column in
#'  `params@species_params`.
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
#' @concept refugeRates
#' @family mizer rate functions
reefEncounter <- function(params, n, n_pp, n_other, t,
                vulnerable = reefVulnerable(params, n, n_pp, n_other, t,
                    new_rd = reefDegrade(params, n, n_pp, n_other, t)),...) {

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

#' Reef feeding level
#'
#' This function replaces the usual [mizerFeedingLevel()] function and returns 
#' the a feeding level of 0 for piscivores.
#'
#' @inheritParams reefEncounter
#' @param encounter A two dimensional array (predator species x predator size) 
#'                  with the encounter rate.
#'
#' @return A two dimensional array (predator species x predator size) with the
#'          feeding level.
#'   
#' @family mizer rate functions
#' @concept extmort
#' @export
reefFeedingLevel <- function(params, n, n_pp, n_other, t, encounter,...) {
    
    # Set predator max intake to Inf
    params@intake_max[params@species_params$satiation == FALSE] <- Inf
    
    # Find mizer feeding level and set to 0 for any NAs (if use wants to set h)
    fl <- mizerFeedingLevel(params, n, n_pp, n_other, t, encounter, ...)
    
    fl[is.na(fl)] <- 0
    
    return(fl)
}

#' Get total predation mortality rate needed to project mizer reef model
#'
#' Calculates the total predation mortality rate \eqn{\mu_{p,i}(w_p)} (in units
#' of 1/year) on each prey species by prey size:
#' 
#' \deqn{\mu_{p.i}(w_p) = \sum_j {\tt pred\_rate}_j(w_p)\, V_{ji}(w_p)\, \theta_{ji}.}{
#'   \mu_{p.i}(w_p) = \sum_j pred_rate_j(w_p) V_{ji}(w_p) \theta_{ji}.}
#'   
#' You would not usually call this function directly but instead 
#' use [getPredMort()], which then calls this function.
#' 
#' @inheritParams reefRates
#' @param vulnerable Array (species x size) with the proportion of individuals 
#'                   that are not protected from predation by refuge
#' @param pred_rate A two dimensional array (predator species x predator size)
#'                  with the predation rate
#'
#' @return A two dimensional array (prey species x prey size) with the 
#'          predation mortality
#' @family mizer rate functions
#' @concept refugeRates
#' @export
reefPredMort <- function(params, n, n_pp, n_other, t, pred_rate,
            vulnerable = reefVulnerable(params, n, n_pp, n_other, t,
            new_rd = reefDegrade(params, n, n_pp, n_other, t)),...) {
    
    # Number of species, number of bins in resource spectrum and full
    # size spectrum
    no_sp <- nrow(params@species_params)
    no_w <- length(params@w)
    no_w_full <- length(params@w_full)
    
    # Get index of species that have grown out of the resource spectrum
    idx_sp <- (length(params@w_full) - 
                   length(params@w) + 1):length(params@w_full)
    pr  <- pred_rate[,idx_sp, drop = FALSE]
    int <- params@interaction
    
    # Find indices of predator species whose foraging is hindered by refuge
    bad_pred  <- which(params@species_params$bad_pred == TRUE)
    good_pred <- which(params@species_params$bad_pred == FALSE)
    
    # Create list of vulnerabilities for each predator
    vul <- vector("list", no_sp)
    vul[bad_pred] <- list(vulnerable)
    vul[good_pred] <- list(matrix(1, nrow = no_sp, ncol = ncol(vulnerable)))
    
    # Loop through predator species to calculate predation mortality on
    # each prey species & size by all predators
    pm <- matrix(0, no_sp, length(params@w))
    dimnames(pm) <- dimnames(vulnerable)
    
    for (i in 1:no_sp){
        i = 1
        # Vulnerability rate of all prey, including resource
        # (species by size) to predator i
        v <- vul[[i]]
        # Predation rate of predator species i on all prey species (by size)
        pr_i <- pr[i,]
        # Predation rate of predator species i (prey species x prey size)
        pr_i <- matrix(rep(pr_i, each = nrow(v)), 
                         nrow = nrow(v), ncol = ncol(v))
        # Interaction of predator species i with all prey
        int_i <- int[i,]
        # Predation rate predator species i * interaction of predator species i
        pr_int <- pr_i * int_i
        # vul*pr_i predation mortality on prey (species by size) by 
        # predator i, adding to pm to sum over all predators
        pm <- pm + v*pr_int
    }

    pred_mort <- pm
    
    return(pred_mort)
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
#' @param n     A matrix of species abundances (species x size).
#' @param n_pp  A vector of the resource abundance by size
#' @param n_other   A list of abundances for other dynamical components of the
#'                  ecosystem
#' @param t     The time for which to do the calculation (Not used by standard
#'              mizer rate functions but useful for extensions with 
#'              time-dependent parameters.)
#' @param ... Unused
#'
#' @return  A named two dimensional array (species x size) with the senescence
#'          mortality rates.
#' @concept extmort
#' @export
reefSenMort <- function(params, n, n_pp, n_other, t, ...) {

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
#' @concept extmort
#' @export
reefMort <- function(params, n, n_pp, n_other, t, f_mort, pred_mort, ...) {
    
    include_sen_mort <- params@other_params$include_sen_mort
    
    if(include_sen_mort == TRUE){
    mizerMort(params, n, n_pp, n_other, t, f_mort, pred_mort, ...) +
        reefSenMort(params, ...)
    } else {
        mizerMort(params, n, n_pp, n_other, t, f_mort, pred_mort, ...)
    }
}




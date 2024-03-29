#' Contribution of unstructured components to the encounter rate
#'
#' The encounter rate \eqn{E_i(w)} for an unstructured resource like algae or
#' detritus is proportional to the total biomass \eqn{B} with a coefficient
#' \eqn{\rho_i(w)} that depends on the functional group \eqn{i} and the size of
#' the consumer:
#'
#' \deqn{E_i(w) = \rho_i(w) B.}
#'
#' The coefficient \eqn{\rho_i(w)} is stored as a matrix (species x size) in
#' the `rho` parameter of the component. It has units 1/year.
#'
#' @param params MizerParams
#' @param n_other Biomasses of unstructured components
#' @param component Name of component whose contribution is requested
#' @param ... Unused
#'
#' @return Array (species x size) with the encounter rate in g/year.
#' @export
#' @concept Uresources
encounter_contribution <- function(params, n_other, component, ...) {
    
    params@other_params[[component]]$rho * n_other[[component]]
}

#' Rescale algae and detritus biomass without changing anything else
#'
#' This multiplies the algae & detritus biomass by a factor and divides the
#' interaction between all groups and unstructured resource by the same
#' factor, so as to keep the total consumption of these resources unchanged.
#' 
#' @param params A MizerParams object
#' @param algae_factor A number to scale algae by
#' @param detritus_factor A number to scale detritus by
#' @return An updated MizerParams object
#' @concept Uresources
#' @export
rescaleComponents <- function(params, algae_factor = 1, detritus_factor = 1) {
    rescale_algae(rescale_detritus(params, detritus_factor),
                  algae_factor)
}

#' Tune unstructured resources with carrying capacities 
#' (algae and detritus) to steady state
#'
#' For models that use unstructured resources with carrying capacities,
#' this functions sets the production rates of detritus and algae so
#' that productions equals consumption at steady state.
#' 
#' With a carrying capacity, the time evolution of the algae biomass 
#' \eqn{B_A(t)} is described by
#'
#'  \deqn{ \frac{dB_A}{dt} = P_A\left( 1 - 
#'                          \frac{B_A}{K_A} \right) - c_A \, B_A }{
#'                 dB_A/dt = P_A * (1 - B_A/ K_A) - c_A * B_A}
#' 
#' where \eqn{K_A} is the system's carrying capacity for algae in grams/ year,
#' \eqn{c_A} is the mass-specific rate of consumption calculated with
#' `algae_consumption()` and \eqn{P_A} is the rate at which algae
#' grows, calculated with `getAlgaeProduction()`.
#' 
#' In this tuning function, the growth of rate of algae is set to 
#' \eqn{(c_A \cdot B_A)/(1-\frac{B_A}{K_A})} grams per meter squared per year
#' so that consumption is equal to production for steady state.
#' 
#' Similarly, the time evolution of the detritus biomass \eqn{B_D(t)} is 
#' described by
#'
#'  \deqn{ \frac{dB_D}{dt} = P_D\left( 1 - 
#'                          \frac{B_D}{K_D} \right) - c_D \, B_D }{
#'                 dB_D/dt = P_D * (1 - B_D/ K_D) - c_D * B_D}
#'                 
#' where \eqn{K_D} is the system's carrying capacity for detritus in grams/ year,
#' \eqn{c_D} is the mass-specific rate of consumption calculated with
#' `detritus_consumption()` and \eqn{P_D} is the rate at which detritus
#' is produced calculated with `getDetritusProduction()`. Total detritus
#' production is given with
#' 
#'  \deqn{p_D = p_{D.f} + p_{D.d} + p_{D.ext}}{
#'        p_D = p_{D.f} + p_{D.d} + p_{D.ext}}
#' 
#' In this tuning function, the external production of detritus is set to 
#' \eqn{(c_D \cdot B_D)/(1-\frac{B_D}{K_D}) - P_{D.f} - P_{D.d}} grams per meter
#' squared per year so that production equals consumption at steady state. 
#'
#' @param params A MizerParams object
#' @param ... unused
#' @return An updated MizerParams object
#' @concept Uresources
#' @seealso [reefSteady()], [algae_dynamics_cc()], [detritus_dynamics_cc()]
#' @export
tuneUR_cc <- function(params,...) {
    
    # algae
    ba <- algae_biomass(params)
    ka <- params@other_params$algae$capacity
    aout <- sum(getAlgaeConsumption(params))
    params@other_params$algae$growth <- (aout*ba)/(1-ba/ka)
    
    # detritus
    params@other_params$detritus$external <- 0
    bd <- detritus_biomass(params)
    kd <- params@other_params$detritus$capacity
    din <- sum(getDetritusProduction(params))
    dout <- sum(getDetritusConsumption(params))
    if (din > dout) {
        warning("The flux of external detritus is negative.")
    }
    params@other_params$detritus$external <- ((dout*bd)/(1-bd/kd)) - din
    
    params
}

#' Tune unstructured resources (algae and detritus) to steady state
#'
#' This first sets the rate of degradation of algae so that for the given
#' abundances, the algae is at steady state. It then sets the rate at which
#' detritus flows in from external sources (e.g. the pelagic zone) so that for
#' the given abundances the detritus is at steady state.
#'
#' @param params A MizerParams object
#' @param ... unused
#' @return An updated MizerParams object
#' @concept Uresources
#' @seealso [reefSteady()], [algae_dynamics()], [detritus_dynamics()]
#' @export
tuneUR <- function(params,...) {

    # algae
    aout <- sum(getAlgaeConsumption(params))
    params@other_params$algae$growth <- aout

    # detritus
    params@other_params$detritus$external <- 0
    din <- sum(getDetritusProduction(params))
    dout <- sum(getDetritusConsumption(params))
    if (din > dout) {
        warning("The flux of external detritus is negative.")
    }
    params@other_params$detritus$external <- dout - din

    params
}

#' Scale reef abundances
#'
#' Multiplies the abundances of all or of selected species by given factors and
#' then retunes the reproductive efficiency accordingly.
#'
#' @details
#' Does not run the system to steady state. For that you should call
#' [reefSteady()] explicitly afterwards.
#'
#' @param params A mizer params object
#' @param factor The factor by which the abundance of each species is multiplied.
#'   This can be specified in two ways:
#'   \itemize{
#'   \item A named numeric vector where the name indicates the species and the
#'     value gives the factor for that species. Only the named species are
#'     affected.
#'   \item  A number that gives the factor for all foreground species.
#'   }
#'
#' @return An object of type \linkS4class{MizerParams}
#' @concept calibration
#' @export
scaleReefAbundance <- function(params, factor) {
    params <- validParams(params)
    assert_that(is.numeric(factor),
                all(factor > 0))
    is_foreground <- !is.na(params@A)
    no_sp <- sum(is_foreground)
    if (length(factor) == 1 && length(names(factor)) == 0) {
        factor <- rep(factor, no_sp)
        names(factor) <- params@species_params$species[is_foreground]
    }
    to_rescale <- names(factor)
    wrong <- setdiff(to_rescale, params@species_params$species)
    if (length(wrong) > 0) {
        stop(paste(wrong, collapse = ", "),
             " do not exist.")
    }
    assert_that(length(to_rescale) == length(factor))
    
    params@initial_n[to_rescale, ] <-
        params@initial_n[to_rescale, ] * factor
    
    return(setBevertonHolt(params, reproduction_level = 1/2))
}

#' Scale model parameters
#'
#' This function scales various model parameters by a given factor.
#'
#' @param params a mizer model object
#' @param factor a numeric value by which to scale the model
#'
#' @return a mizer model object with scaled parameters
#' @concept calibration
#' @export
scaleReefModel <- function(params, factor) {

    # Algae
    params@other_params[["algae"]]$rho <-
        params@other_params[["algae"]]$rho / factor
    params@species_params$rho_algae <-
        params@species_params$rho_algae / factor
    params@other_params$algae$growth <-
        params@other_params$algae$growth * factor

    # Detritus
    params@other_params[["detritus"]]$rho <-
        params@other_params[["detritus"]]$rho / factor
    params@species_params$rho_detritus <-
        params@species_params$rho_detritus / factor
    params@other_params$detritus$external <-
        params@other_params$detritus$external * factor
    
    # now comes the code of mizer's standard scaleModel()
    params <- validParams(params)
    assert_that(is.number(factor), factor > 0)
    params@cc_pp <- params@cc_pp * factor
    params@resource_params$kappa <- params@resource_params$kappa * factor
    if ("r_max" %in% names(params@species_params)) {
        params@species_params$R_max <- params@species_params$r_max
        params@species_params$r_max <- NULL
        message("The 'r_max' column has been renamed to 'R_max'.")
    }
    if ("R_max" %in% names(params@species_params)) {
        params@species_params$R_max <- params@species_params$R_max * factor
    }
    params@search_vol <- params@search_vol/factor
    if ("gamma" %in% names(params@species_params)) {
        params@species_params$gamma <- params@species_params$gamma/factor
    }
    initial_n_other <- params@initial_n_other
    for (res in names(initial_n_other)) {
        initial_n_other[[res]] <- initial_n_other[[res]] * factor
    }
    initialN(params) <- params@initial_n * factor
    initialNResource(params) <- params@initial_n_pp * factor
    initialNOther(params) <- initial_n_other
    params@sc <- params@sc * factor
    return(params)
}

#' Scale background down by a factor
#' 
#' Replaces scale down background function
#' 
#' @param params a mizer model object
#' @param factor A number giving the factor by which the background abundance
#'   will be reduced
#' @export
#' @concept Uresources
scaleReefBackground <- function(params, factor) {
    scaleReefAbundance(params, factor = factor) %>%
        scaleModel(factor = 1 / factor)
}



#' Calibrate the scale of a mizerReef model to match total observed biomass
#' 
#' This function replaces mizer's calibrateBiomass function to include 
#' unstructured resources. Given a Mizer Params object, it returns an updated 
#' MizerParams object which is rescaled with [scaleReefModel()] so that the 
#' total biomass in the model agrees with the total observed biomass.
#' 
#' Biomass observations usually only include individuals above a certain size.
#' This size should be specified in a biomass_cutoff column of the species
#' parameter data frame. If this is missing, it is assumed that all sizes are
#' included in the observed biomass, i.e., it includes larval biomass.
#' 
#' After using this function the total biomass in the model will match the
#' total biomass, summed over all species. However the biomasses of the
#' individual species will not match observations yet, with some species
#' having biomasses that are too high and others too low. So after this
#' function use the mizer function matchReefBiomasses() to match the 
#' biomasses for each group. 
#' 
#' @param params A MizerParams object
#' @return A MizerParams object
#' @concept calibration
#' @export
calibrateReefBiomass <- function(params) {
    if ((!("biomass_observed" %in% names(params@species_params))) ||
        all(is.na(params@species_params$biomass_observed))) {
        return(params)
    }
    no_sp <- nrow(params@species_params)
    cutoff <- params@species_params$biomass_cutoff
    # When no cutoff known, set it to 0
    if (is.null(cutoff)) cutoff <- rep(0, no_sp)
    cutoff[is.na(cutoff)] <- 0
    observed <- params@species_params$biomass_observed
    observed_total <- sum(observed, na.rm = TRUE)
    sp_observed <- which(!is.na(observed))
    model_total <- 0
    for (sp_idx in sp_observed) {
        model_total <- 
            model_total + 
            sum((params@initial_n[sp_idx, ] * params@w * params@dw)
                [params@w >= cutoff[[sp_idx]]])
    }
    scaleReefModel(params, factor = observed_total / model_total)
}

# The following is a copy of the code for `calibrateBiomass()` just with
# the text replacements "Biomass" -> "Number" and "biomass" to "number" and
# the removal of the `params@w` factor in the calculations.

#' Calibrate the model scale to match total observed number
#'
#' Replaces mizer's [mizer::calibrateNumber()] function. Given a MizerParams
#' object `params` for which number observations are available for at least 
#' some species via the `number_observed` column in the species_params data 
#' frame, this function returns an updated MizerParams object which is 
#' rescaled with [scaleReefModel()] so that the total number in
#' the model agrees with the total observed number.
#'
#' Number observations usually only include individuals above a certain size.
#' This size should be specified in a number_cutoff column of the species
#' parameter data frame. If this is missing, it is assumed that all sizes are
#' included in the observed number, i.e., it includes larval number.
#'
#' After using this function the total number in the model will match the
#' total number, summed over all species. However the numbers of the
#' individual species will not match observations yet, with some species
#' having numbers that are too high and others too low. So after this
#' function you may want to use [matchNumbers()]. This is described in the
#' blog post at https://bit.ly/2YqXESV.
#'
#' If you have observations of the yearly yield instead of numbers, you can
#' use [calibrateYield()] instead of this function.
#'
#' @param params A MizerParams object
#' @return A MizerParams object
#' @concept calibration
#' @export
calibrateReefNumber <- function(params) {
    if ((!("number_observed" %in% names(params@species_params))) ||
        all(is.na(params@species_params$number_observed))) {
        return(params)
    }
    no_sp <- nrow(params@species_params)
    cutoff <- params@species_params$number_cutoff
    # When no cutoff known, set it to 0
    if (is.null(cutoff)) cutoff <- rep(0, no_sp)
    cutoff[is.na(cutoff)] <- 0
    observed <- params@species_params$number_observed
    observed_total <- sum(observed, na.rm = TRUE)
    sp_observed <- which(!is.na(observed))
    model_total <- 0
    for (sp_idx in sp_observed) {
        model_total <-
            model_total +
            sum((params@initial_n[sp_idx, ] * params@dw)
                [params@w >= cutoff[[sp_idx]]])
    }
    scaleReefModel(params, factor = observed_total / model_total)
}


#' Hold resource dynamics constant
#' 
#' @param params MizerParams object
#' @param n_other Biomasses of unstructured components
#' @param component Name of component to view dynamics for
#' @param ... Unused
#' @concept Uresources
#' @export
constant_dynamics <- function(params, n_other, component, ...) {
    n_other[[component]]
}

#' Match observed growth rates
#' 
#' This does the same as `mizer::matchGrowth()` but in addition also rescales
#' the consumption rates of algae and detritus.
#' 
#' @param params A MizerParams object
#' @param species The species to be affected. Optional. By default all species
#'   for which growth information is available will be affected. A vector of
#'   species names, or a numeric vector with the species indices, or a logical
#'   vector indicating for each species whether it is to be affected (TRUE) or
#'   not.
#' @param keep A string determining which quantity is to be kept constant. The
#'   choices are "egg" which keeps the egg density constant, "biomass" which 
#'   keeps the total biomass of the species constant and "number" which keeps
#'   the total number of individuals constant.
#' @return A modified MizerParams object with rescaled rates and rescaled 
#'   species parameters `gamma`,`h`, `ks` and `k`.
#' @concept calibration
#' @export
matchReefGrowth <- function(params, species = NULL,
                            keep = c("egg", "biomass", "number")) {
    
    assert_that(is(params, "MizerParams"))
    sel <- valid_species_arg(params, species = species, 
                             return.logical = TRUE)
    sp <- params@species_params
    keep <- match.arg(keep)
    
    biomass <- getBiomass(params)
    number <- getN(params)
    
    sp <- set_species_param_default(sp, "age_mat", NA)
    # If age at maturity is not specified, calculate it from von Bertalanffy
    if (all(c("k_vb", "w_inf") %in% names(sp))) {
        age_mat_vB <- age_mat_vB(params)
        sp <- mizer::set_species_param_default(sp, "age_mat", age_mat_vB)
    }
    
    # Don't affect species where no age at maturity is available
    sel <- sel & !is.na(sp$age_mat)
    
    factor <- age_mat(params)[sel] / sp$age_mat[sel]
    
    params@search_vol[sel, ] <- params@search_vol[sel, ] * factor
    params@intake_max[sel, ] <- params@intake_max[sel, ] * factor
    params@metab[sel, ] <- params@metab[sel, ] * factor
    params@species_params$gamma[sel] <- sp$gamma[sel] * factor
    params@species_params[sel, "h"] <- sp[sel, "h"] * factor
    if ("ks" %in% names(sp)) {
        params@species_params$ks[sel] <- sp$ks[sel] * factor
    }
    if ("k" %in% names(sp)) {
        params@species_params[sel, "k"] <- sp[sel, "k"] * factor
    }
    
    # rescale consumption of algae and detritus
    
    params@other_params$algae$rho[sel, ] <- 
        params@other_params$algae$rho[sel, ] * factor
    params@species_params$rho_algae[sel] <-
        params@species_params$rho_algae[sel] * factor
    
    params@other_params$detritus$rho[sel, ] <- 
        params@other_params$detritus$rho[sel, ] * factor
    params@species_params$rho_detritus[sel] <-
        params@species_params$rho_detritus[sel] * factor
    
    params <- steadySingleSpecies(params, species = sel)
    
    if (keep == "biomass") {
        factor <- biomass / getBiomass(params)
        params@initial_n <- params@initial_n * factor
    }
    if (keep == "number") {
        factor <- number / getN(params)
        params@initial_n <- params@initial_n * factor
    }
    
    setBevertonHolt(params)
}

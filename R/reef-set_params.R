#' Checks unstructured resource interaction matrix
#'
#' @section Adding unstructured resources:
#'
#'      mizerReef supports two resource spectra that are not size- structured.
#'      Algae is consumed by herbivorous fish, while detritus is consumed by
#'      herbivorous fish and benthic invertebrates.
#'
#'      The resource interaction matrix \eqn{\theta_{ki} modifies the
#'      interaction of each functional group \eqn{i} with each resource \eqn{k}
#'      in the model. This can be used for example to allow for different
#'      spatial overlap among the species. This function checks if provided
#'      unstructured resource interaction matrices are valid.
#'
#' @param params MizerParams object
#' @param UR_interaction Interaction matrix for unstructured resources
#'  (species x resource)
#'
#' Optional parameters:
#' @param algae_growth The initial growth rate of algae in grams/year/m^-2.
#'  This value is reset to match consumption in the `[reefSteady()]`  function
#'  so that steady state abundances match given values.
#' @param prop_decomp The proportion of waste material that decomposes to
#'  become part of the detritus pool.
#' @param d.external The rate at which detritus biomass sinks from the pelagic
#'  zone and becomes part of the detritus pool in grams per year.
#'
#' @return `setUResourceParams`: MizerParams object with updated unstructured
#'  resource parameters
#' @export
#' @family functions for setting parameters
setURParams <- function(params,
                        UR_interaction,
                        algae_growth = NULL,
                        prop_decomp = NULL,
                        d.external = NULL) {

    # Check if mizerParams is valid
    assert_that(is(params, "MizerParams"))

    # Find number of species for checks
    no_sp = nrow(params@species_params)

    # check if matrix has right names
    res_cols <- c('algae','detritus')
    if(!setequal(names(UR_interaction), res_cols)){
        stop("Uresource_interaction needs to have columns named 'algae' and
             'detritus'.")
    }

    # check if detritus and algae interactions are right length
    if(length(UR_interaction$algae) != no_sp) {
        stop("The 'algae' column should have a value for every functional
             group.")
    }
    if(length(UR_interaction$detritus) != no_sp){
        stop("The 'detritus' column should have a value for every functional
             group.")
    }

    # Check if values are numeric
    if(!all(sapply(UR_interaction, is.numeric))) {
        stop("The ntries for algae & detritus interaction should be numeric.")
    }

    # Check if values are between 0 and 1
    if(!is.matrix(UR_interaction)) {
        UR <- as.matrix(UR_interaction)

        if(any(UR < 0)) {
            stop("Entries for algae & detritus interaction should be
                 non-negative.")
        }
        if(any(UR > 1 )) {
            stop("Entries for algae & detritus interaction should be
                 between 0 and 1.")
        }
    }

    # Set default algae growth rate
    if(is.null(algae_growth)){ params@other_params$algae_growth <- 200
    } else{ params@other_params$algae_growth <- algae_growth }

    # Set default proportion of waste that becomes part of the detritus pool
    if(is.null(prop_decomp)){ params@other_params$prop_decomp <- 0.2
    } else { params@other_params$prop_decomp <- prop_decomp }

    # Set default external detritus
    if(is.null(d.external)){ params@other_params$d.external <- 0.1
    } else { params@other_params$d.external <- d.external }

    # Add values as columns to species params data frame
    params@species_params$interaction_algae <- UR_interaction$algae
    params@species_params$interaction_detritus <- UR_interaction$detritus

    # Note time sim was modified
    params@time_modified <- lubridate::now()

    return(params)
}

#' Calculate rho parameter values using a critical feeding level
#'
#' This function uses the `Uresource` interaction matrix as well as a
#' given critical feeding level to find \eqn{\rho_i} for all unstructured
#' resources. It then stores these values in the `params@species_params`
#' data frame. There are two unstructured resources in a mizerReef model,
#' algae and detritus.
#'
#' @section The \eqn{\rho_{ki}} Parameter:
#'
#'      For each consumer group \eqn{i} and unstructured resource \eqn{k},
#'      a parameter \eqn{\rho_{ki}} determines the rate at which individuals of
#'      that functional group encounter unstructured resources. The parameters
#'      \eqn{\rho_{ki}} have units of \eqn{g^{-n}} per year. They are non-zero
#'      only for species that forage on resource \eqn{k}.
#'
#'      The rate is assumed to scale with the size of the predator raised to
#'      an allometric exponent \eqn{n} which is taken to be the same as the
#'      scaling exponent of the maximum intake rate for consumers. The
#'      encounter rate for unstructured resources is given by:
#'
#'      \deqn{ E_{i}(w)= \rho_{ki} w^n B_{k} }{ \rho_{ki} w^n B_{k} }
#'
#'      where \eqn{B_{k}} is the biomass of unstructured resource \eqn{k}.
#'
#' @section Using critical feeding level to determine \eqn{\rho_{ki}}:
#'
#'      This function finds reasonable values for \rho_{ki} parameters given
#'      a critical feeding level.
#'
#'      First, it set the critical feeding level  to the provided value
#'      for fish at the maximum size. Then it finds the rate at which each the
#'      functional group is currently encountering food.
#'
#'      \eqn{\rho_{ki}} is then set so that if the current encounter rate is
#'      too low to meet satiation (critical feeding level \eqn{f0}), the group
#'      feeds on resource \eqn{k} to make up the difference.
#'
#' @param params MizerParams object
#' @param UR_interaction Interaction matrix for unstructured resources
#'  (species x resource)
#'
#' Optional parameters:
#' @param crit_feed the desired critical feeding level for herbivores and
#'  detritivores
#'
#' @return `setUResourceParams`: MizerParams object with updated unstructured
#'  resource parameters
#' @export
#' @family functions for setting parameters
setRho <- function(params,
                   crit_feed = NULL,
                   n = 3/4, ...){
    # Determine the necessary detritus and algae encounter rates so that at
    #maximum size the group has feeding level f0
    if(is.null(crit_feed)){ crit_feed <- 0.6 }
    f0 <- set_species_param_default(params@species_params, "f0", crit_feed)$f0

    # Get interaction of each species with detritus and algae
    ia <- params@species_params$interaction_algae
    id <- params@species_params$interaction_detritus

    # Calculate encounter rates divided by w^n of largest individuals
    E <- getEncounter(params)[, length(params@w)] /
        (params@w[length(params@w)] ^ n)

    # Calculate rho for each unstructured resource
    # f0*h/(1-f0) is the encounter rate when feeding level is f0
    # We subtract E so that if feeding level is too low, they eat
    # algae to replace it. For unstructured resources
    # encounter rate = rho * w^n * B_A, and multiply by interaction
    rho_alg <- pmax(0, f0 * params@species_params$h / (1 - f0) - E) * ia
    rho_det <- pmax(0, f0 * params@species_params$h / (1 - f0) - E) * id

    # Store new rho values in species_params data frame
    params@species_params$rho_algae    <- rho_alg
    params@species_params$rho_detritus <- rho_det

    # Note time sim was modified
    params@time_modified <- lubridate::now()

    return(params)
}


#' Set the parameters for external mortality
#'
#' Checks if given parameters are valid, sets defaults, and then stores them in
#' the mizerParams object.
#'
#' mizerReef models contain two sources of external mortality, external
#' predation and senescence.
#'
#' @section External predation:
#'
#' The external predation accounts for all the mortality that is not due to
#' predation by predators included in the model. It is a rate with units
#'  1/year.
#'
#' External predation mortality is assumed to decrease allometrically with
#' body size, given by:
#'
#' \eqn{\mu_{ext}(w) = z_{0.i}}. The value of the constant \eqn{z_0} for each
#' species is taken from the `z0` column of the species parameter data frame, if
#' that column exists. Otherwise it is calculated as
#' \deqn{\mu_{ext}(w) = {\tt z0pre}_i\, w_{inf}^{\tt z0exp}.}{z_{0.i} = z0pre_i w_{inf}^{z0exp}.}
#'
#' @section Senescence mortality:
#'
#'      Senescence mortality is implemented by default in mizerReef. You do not
#'      need to set these parameters. This function should only be used to
#'      change default values.
#'
#'      Senescence mortality \eqn{\mu_{sen.i}(w)} is used to represent
#'      mortality caused by external sources such as illness or age. This is
#'      addition to external mortality, \eqn{\mu_{ext.i}(w)}, which represents
#'      all mortality that is not due to fishing or predation by predators
#'      included in the model. The rate of senescence mortality is given by:
#'
#'      \deqn{\mu_{sen.i}(w) = k_{sen}\left(\frac{w}{w_{max.i}}\right)^{p_{sen}}}
#'           {\mu_{sen.i}(w) = k_{sen}(w/w_{max.i})^{p_{sen}}}
#'
#'      where \eqn{k_{sen}} is the rate of senescence mortality, \eqn{p_sen}
#'      defines the slope of the senescence curve, \eqn{w_max.i} maximum body
#'      size of group \eqn{i} in grams.
#'
#' @param params MizerParams object
#'
#' Optional parameters:
#' @param ext_mort_params Named list containing desired mortality parameters
#' @return `setExtMortParams`: MizerParams object with updated mortality
#'  senescence parameters
#' @export
#' @family functions for setting parameters
setExtMortParams <- function(params,
                             ext_mort_params = NULL) {

    # Check if mizerParams is valid
    assert_that(is(params, "MizerParams"))

    # Check if user provided valid mortality parameters
    if(!is.null(ext_mort_params)) {
        if (!all(sapply(ext_mort_params, is.numeric))) {
            stop("The external mortality parameters should be numeric.")
        }
        if(!is.matrix(ext_mort_params)) {
            ext_mort_params <- as.matrix(ext_mort_params)
            if(!all(ext_mort_params >= 0 )){
                stop("The external mortality parameters should be
                     nonnegative")
            }
        }
        if(!("nat_mort" %in% colnames(ext_mort_params))){
            stop("The external mortality parameters dataframe needs a column
                 called 'nat_mort' with the residual natural mortality rate.")
        }
        if(!("sen_prop" %in% colnames(ext_mort_params))){
            stop("The external mortality parameters dataframe needs a column
                 called 'sen_prop' with the external mortality rate.")
        }
        if(!("sen_curve" %in% colnames(ext_mort_params))){
            stop("The external mortality parameters dataframe needs a column
                 called 'sen_curve' with the exponent for the external
                 curve.")
        }
    } else {

        ext_mort_params <- vector('list',3)
        names(ext_mort_params) <- c('nat_mort','sen_prop','sen_curve')

        # sen_mort_params <- vector('list',3)
        # names(sen_mort_params) <- c('sen_prop','sen_curve','sen_length')
        # External predation
        ext_mort_params$nat_mort  <- 0.2
        # senescence
        ext_mort_params$sen_prop  <- 0.1
        ext_mort_params$sen_curve <- 0.3
        #sen_mort_params$sen_length <- 50
    }

    # Store in params
    params@other_params[['ext_mort_params']] <- ext_mort_params

    params@time_modified <- lubridate::now()

    return(params)
}

#' Set the refuge profile parameters
#'
#' Checks if given refuge parameters are valid and then stores them in the
#' mizerParams object
#'
#' @section Setting the refuge profile:
#'
#'  Refuge profiles account for the protective behavior of prey living in
#'  high complexity environments (e.g. coral reefs) with access to predation
#'  refuge. The refuge profile defines the proportion of fish within
#'  user-defined length bins that are protected from being encountered
#'  by a predator.
#'
#'  A unique refuge profile is generated for each predator x prey combination
#'  based on the given refuge profile parameters as well as four values from
#'  `params@species_params`: length to weight conversion values `a` and `b`,
#'  `refuge_user`, which is true for groups utilize that predation refuge,
#'  and `bad_pred`, which is true for predators whose foraging is hindered
#'  by refuge.
#'
#'  The refuge profile is used when calculating the food encounter rate in
#'  [reefEncounter()] and the predation mortality rate in [reefPredMort()]. Its
#'  entries are dimensionless numbers. If no refuge is available then
#'  predator-prey interactions are determined entirely by size-preference.
#'
#'  The mizerReef package provides three methods to define the refuge profile:
#'
#'      1. ``simple`` -  This method is preferred for data-poor reefs or reefs
#'      where the refuge distribution is unknown. It is also ideal for systems
#'      where only one functional group is expected to be utilizing refuge. The
#'      proportion of fish with access to refuge \eqn{ R_j(w_p) }$ is given by:
#'
#'      \deqn{ R_j(w_p) =
#'          \frac{-ref}{1 + e^{\left(-\alpha (w - W_{max})\right)}} + ref }
#'          {R_j(w_p) = -r/(1 + e^{(-\alpha (w - W_{max}))} + ref }
#'
#'      where \eqn{W_{max}} defines the prey body size at which no crevices in
#'      the reef are large enough to act as a refuge. Refuge is available to a
#'      constant proportion \eqn{ref} of fish smaller than \eqn{W_{max}}. The
#'      slope \eqn{\alpha} describes the sharpness of the cutoff for fish
#'      larger than \eqn{W_{max}}. The default value for \eqn{\alpha} sets a
#'      steep slope of \eqn{100}.
#'
#'      For this method, `method_params` should contain columns named
#'      `prop_protect` and `max_L` that give the proportion of fish to protect
#'      and the maximum length protected (cm), respectively.
#'
#'      2. ``binned`` - This method is appropriate for theoretical applications
#'      and does not rely on empirical data. It sets refuge to a constant
#'      proportion of fish within a given size range. The  proportion of fish
#'      in group \eqn{j} with access to refuge is given by:
#'
#'      \deqn{ R_j(w_p) = r_k ~~~~~~~ w_p ∈ (~w_{k-1}, w_k~] }
#'      {R_j(w_p) = r_k ~~~~~~~ w_p ∈ (~w_{k-1}, w_k~]
#'
#'      where \eqn{r_k} is the proportion of fish with access to refuge in
#'      size class \eqn{k}.
#'
#'      For this method, `method_params` should contain columns named
#'      `start_L`and `end_L` which contain the starting and ending lengths (cm)
#'      of each size bin and `prop_protect`, the proportion protected within
#'      each corresponding size bin.
#'
#'      3. ``data`` -  This method is appropriate when data on the number of
#'      refuge holes present within defined fish length bins is available. The
#'      proportion of fish in size class \eqn{k} with access to refuge is given
#'      by:
#'
#'      \deqn{R_j(w_p) = \tau\cdot\frac{\eta_k}{\sum_{i}\int_{w_{k-1}}^{w_k}
#'      N_i(w)~dw}  ~~~~~~~ w_p ∈ (~w_{k-1}, w_k~]}
#'      {R_j(w_p) = \tau \eta_k/(\sum_{i}\int_{w_{k-1}}^{w_k}
#'      N_i(w)~dw)  ~~~~~~~ w_p ∈ (~w_{k-1}, w_k~]}
#'
#'      where \eqn{\tau} is the proportion of fish with access to refuge that
#'      are expected to actually utilize  it, \eqn{\eta_k} is the density of
#'      refuges in size range \eqn{(w_{k-1}, w_k]} and
#'       \eqn{\sum_{i}\int_{w_{k-1}}^{w_k} N_i(w)~dw} gives the total
#'      number of fish from any group in size range \eqn{(w_{k-1}, w_k]}.
#'      This represents the density of competitors for refuges in
#'      size class \eqn{k}.
#'
#'      For this method, `method_params` should contain columns named
#'      `start_L`and `end_L` which contain the starting and ending lengths (cm)
#'      of each size bin and `refuge_density`, the number of refuges available
#'      in each size bin.
#'
#'  This function checks that the supplied refuge parameters are valid, adds
#'  relevant columns to the `species_params` dataframe, and stores refuge
#'  parameters in the `other_params` slot of the `params` object.
#'
#'  Refuge profile parameters can be input in a spreadsheet program and saved
#'  as a .csv file. The data can then be read into R using the command
#'  `read.csv()`.
#'
#' @param params MizerParams object
#' @param method The desired method for setting up benthic refuge,
#'  can be "simple", "binned", or "data"
#' @param method_params A data frame containing values specific to each
#' method for calculating refuge
#' @param refuge_user A vector of logical values indicating whether each
#' functional group uses refuge, TRUE indicates the group uses refuge while
#' false indicates that they do not. Alternatively,t hese can be included in
#' the `params@species_params` dataframe.
#' @param bad_predator A vector of logical values indicating whether hunting is
#' inhibited by refuge for this functional group. FALSE indicates that this
#' species is able to encounter prey within refuge (e.g. eels). Alternatively,
#' these can also be included in the `params@species_params` dataframe.
#' @param pisc A vector of logical values indicating whether this
#' functional group are generally piscivores. FALSE indicates that this group
#' feeds only from the resource spectra. Alternatively,
#' these can also be included in the `params@species_params` dataframe.
#'
#' Optional parameters:
#' @param w_settle  The body weight (g) at which fish settle onto the reef.
#'                  Fish smaller than this are considered to be larval and
#'                  thus too small to use predation refuge
#' @param max_protect The maximum proportion of fish (any size class)
#'                      protected by refuge
#' @param tau The proportion of fish with access to refuge that actually use it
#'
#' @return `setRefuge`: A MizerParams object with updated refuge parameters
#' @export
#' @family functions for setting parameters
setRefuge <- function(params,
                      method,
                      w_settle = NULL,
                      max_protect = NULL,
                      tau = NULL,
                      method_params,
                      refuge_user = NULL,
                      bad_predator = NULL,
                      pisc = NULL,...) {

    # Check if mizerParams is valid
    assert_that(is(params, "MizerParams"))

    # Check that a and b parameters are present for all species -
    # needed for l2w conversion
    if (anyNA(params@species_params[["a"]]) ||
        anyNA(params@species_params[["b"]])) {
        stop("There must be no NAs in the species_params columns 'a' and 'b'.")
    }

    # Find number of species for checks
    no_sp = nrow(params@species_params)

    # Check if the user provided one of the available methods
    method_options <- c('simple','binned','data')
    if(is.null(method)) {
        stop("You must provide the method to calculate the refuge profile.")
    } else if(!is.element(method, method_options)) {
        stop("Method must be 'simple','binned', or 'data'.")
    }

    # Set default values for parameters used by all methods

    # Minimum size of fish protected by refuges at measured scale
    if(is.null(w_settle)){
        w_settle <- 0.01
        # # Calculate minimum weight that can use refuge for each group
        # w_settle <- params@species_params[["a"]] *
        #     min_ref_length ^ params@species_params[["b"]]
    } else {
        if(!is.numeric(w_settle)) {
            stop("w_settle should be numeric.")
        }
        if(w_settle < 0) {
        stop("w_settle must be non-negative.")
        }
    }

    # Maximum proportion of fish protected by refuge
    if(is.null(max_protect)){
        max_protect <- 0.98
    } else {
        if(!is.numeric(max_protect)) { stop("max_protect should be numeric.") }
        if(max_protect < 0 || max_protect > 1) {
            stop("max_protect should be a proportion between 0 and 1")
        }
    }

    # Proportion of fish with access to refuge that are expected to utilize it
    if(is.null(tau)){
        tau <- 0.90
    } else {
        if(!is.numeric(tau)) { stop("tau should be numeric.") }
        if(refuge_params$tau < 0 || refuge_params$tau > 1) {
            stop("tau should be a proportion between 0 and 1")
        }
    }

    # Store all values in refuge_params data frame
    refuge_params <- data.frame(method, w_settle, max_protect, tau)

    # check if method_params values are positive and numeric
    if (!is.matrix(method_params)) {
        rmp <- as.matrix(method_params)
        # Check that values are numbers
        if (!all(sapply(rmp, is.numeric))) {
            stop("The method parameters should be numeric.")
        }
        # Check that all values of refuge method matrix are positive
        if (!all(rmp >= 0)) {
            stop("The method parameters must be non-negative.")
        }
    }

    cnames = colnames(method_params)

    # Check names of method_params for simple method
    if (refuge_params$method == "simple") {
        if(!("prop_protect" %in% cnames)) {
            stop("The simple method parameters dataframe needs a column called
                 'prop_protect' with the proportion of fish protected.")
        } else if(method_params$prop_protect < 0 ||
                  method_params$prop_protect > 1) {
            stop("prop_protect should be a proportion between 0 and 1")
        }
        if(!("max_L" %in% cnames)) {
            stop("The simple method parameters dataframe needs a column called
                 'max_L' with the maximum length (cm) of protected fish.")
        }
        if (!is.null(method_params$slope)){ slope <- 100 }
    }

    # Check names of method_params for binned method
    if (refuge_params$method == "binned") {
        if(!("start_L" %in% cnames)) {
            stop("The binned method parameters dataframe needs a column called
                 'start_L' with the starting lengths (cm) for each size bin.")
        }
        if(!("end_L" %in% cnames)) {
            stop("The binned method parameters dataframe needs a column called
                 'end_L' with the end lengths (cm) for each size bin.")
        }
        if(!("prop_protect" %in% cnames)) {
            stop("The binned method parameters dataframe needs a column called
                 'prop_protect' with the proportion of fish protected
                 for each length bin.")
        } else if(any(method_params$prop_protect < 0) ||
                  any(method_params$prop_protect > 1)) {
            stop("prop_protect should be a proportion between 0 and 1")
        }
        if (!all(method_params$start_L < method_params$end_L)) {
            stop("All bin start lengths must be less than bin end lengths.")
        }
    }

    # Check names of method_params for binned method
    if (refuge_params$method == "data") {
        if(!("start_L" %in% cnames)) {
            stop("The data method parameters dataframe needs a column called
                 'start_L' with the starting lengths (cm) for each size bin.")
        }
        if(!("end_L" %in% cnames)) {
            stop("The data method parameters dataframe needs a column called
                 'end_L' with the end lengths (cm) for each size bin.")
        }
        if(!("refuge_density" %in% cnames)) {
            stop("The data method parameters dataframe needs a column called
                 'refuge_density' with the proportion of fish protected
                 for each bin.")
        }
        if (!all(method_params$start_L < method_params$end_L)) {
            stop("All bin start lengths must be less than bin end lengths.")
        }
    }

    # Check that refuge_user, bad_pred, and pisc are present,
    # the right length, and logical values
    if(!('refuge_user' %in% colnames(params@species_params))){
        if(is.null(refuge_user)){
            stop("You need to provide values for refuge_user")
        } else if (!is.logical(refuge_user)) {
            stop("The refuge_user values should be logical.")
        }
        if(length(refuge_user) != no_sp) {
        stop("refuge_user should have a value for every group.")
        }
        params@species_params$refuge_user <- refuge_user
    }

    if(!('bad_pred' %in% colnames(params@species_params))){
        if(is.null(bad_pred)){
            stop("You need to provide values for bad_pred")
        } else if (!is.logical(bad_pred)) {
            stop("The bad_predr values should be logical.")
        }
        if(length(refuge_user) != no_sp) {
            stop("bad_pred should have a value for every group.")
        }
        params@species_params$bad_pred <- bad_pred
    }

    if(!('pisc' %in% colnames(params@species_params))){
        if(is.null(pisc)){
            stop("You need to provide values for pisc")
        } else if (!is.logical(pisc)) {
            stop("The pisc values should be logical.")
        }
        if(length(refuge_user) != no_sp) {
            stop("bad_pred should have a value for every group.")
        }
        params@species_params$pisc <- pisc
    }

    # Store in params
    params@other_params[['refuge_params']] <- as.data.frame(refuge_params)

    params@other_params[['method_params']] <- as.data.frame(method_params)

    params@species_params$w_settle <- w_settle

    params@time_modified <- lubridate::now()

    return(params)
}

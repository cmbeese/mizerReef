#' Checks unstructured resource interaction matrix
#'
#' @section Adding unstructured resources:
#'
#'      mizerReef supports two resource spectra that are not size- structured.
#'      Algae is consumed by herbivorous fish, while detritus is consumed by
#'      herbivorous fish and benthic invertebrates. This function sets the 
#'      interaction matrix for these resources as well as any default
#'      parameters necessary to structure them.
#'
#'      The resource interaction matrix \eqn{\theta_{ki} modifies the
#'      interaction of each functional group \eqn{i} with each unstructured
#'      resource \eqn{k} in the model. This can be used for example to allow 
#'      for different diet preferences on each unstructured resource. 
#'      
#'      Note that interaction with size structured resources, such as
#'      plankton, is still set with the `resource_interaction` column of
#'      the species parameters dataframe. 
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
#'  zone and becomes part of the detritus pool in grams per year. This value 
#'  is reset to make up any differences in consumption and production in 
#'  the `[reefSteady()]` function so that steady state abundances match 
#'  observed values.
#'
#' @return `setUResourceParams`: MizerParams object with updated unstructured
#'  resource parameters
#' @export
#' @concept Unstructured resources
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
        stop("UR_interaction needs to have columns named 'algae' and
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
        stop("The entries for algae & detritus interaction should be numeric.")
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
    if(is.null(algae_growth)){ params@other_params$algae_growth <- 1000
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

#' Set the parameters for external mortality
#'
#' mizerReef models contain two sources of external mortality, residual
#' natural mortality and senescence. This function checks if given parameters 
#' are valid, sets defaults, and then stores them in the mizerParams object.
#' 
#' External mortality implemented by default in mizerReef. You do not need to 
#' set these parameters. This function should only be used to change default
#' values.
#'
#' @section Residual natural mortality:
#' 
#'      Residual natural mortality accounts for any external predation or
#'      fishing mortality that is not explicitly included in the model. It is
#'      assumed to decrease allometrically with body size. Residual natural
#'      mortality is a rate with units 1/year given by:
#'      
#'      \deqn{\mu_{nat.i}(w) = \mu_{nat}\, w^{1-n}.}
#'           {\mu_{nat.i}(w) = \mu_{nat}\, w^{1-n}.}
#'           
#'       Here \eqn{\mu_{nat}} is the residual natural mortality rate at size
#'       1 g and \eqn{n} is the allometric scaling exponent. In mizerReef, 
#'       these default to \eqn{\mu_{nat} = 0.2} and \eqn{n = 0.75}.
#'
#' @section Senescence mortality:
#'
#'      Senescence mortality \eqn{\mu_{sen.i}(w)} is used to represent
#'      mortality caused by background sources such as illness or age. The 
#'      rate of senescence mortality (in 1/year) is given by:
#'
#'      \deqn{\mu_{sen.i}(w) = k_{sen}\left(
#'                              \frac{log_{10}(w)}{log_{10}(w_{max.i})}
#'                              \right)^{p_{sen}}}
#'           {\mu_{sen.i}(w) = k_{sen}
#'           (log_{10}(w)/log_{10}(w_{max.i}))^{p_{sen}}}
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
#' @concept External mortality
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
#'  A unique refuge profile is generated for each predator group x 
#'  prey group x prey size combination based on the given refuge profile
#'  parameters as well as four values from `params@species_params`: length 
#'  to weight conversion values `a` and `b`, `refuge_user`, which is true 
#'  for groups utilize that predation refuge, and `bad_pred`, which is false 
#'  for predator groups whose body shape or predatory strategy allow them to
#'  access fish within refuge (e.g. eels).
#'  
#'  To ensure some food is always available to predators, the maximum 
#'  proportion of fish protected by refuge in any size class is set by
#'  `max_protect`.
#'
#'  The refuge profile is used when calculating the food encounter rate in
#'  [reefEncounter()] and the predation mortality rate in [reefPredMort()]. 
#'  Its entries are dimensionless values between 0 and 1 which represent the
#'  proportion of fish in the corresponding prey and size categories that are
#'  hidden within refuge and thus cannot be encountered by predators. If no
#'  refuge is available then predator-prey interactions are determined 
#'  entirely by size-preference.
#'
#'  The mizerReef package provides three methods to define the refuge profile:
#'
#'      1.``sigmoidal`` - This method is preferred for data-poor reefs or reefs
#'      where the refuge distribution is unknown. It is also ideal for systems
#'      where only one functional group is expected to be utilizing refuge. The
#'      proportion of fish with access to refuge \eqn{ R_j(w_p) }$ is given by:
#'
#'      \deqn{ R_j(w_p) =
#'          \frac{r}{1 + e^{\left(\alpha (w - W_{refuge})\right)}}}
#'          {R_j(w_p) = r_{peak}/(1 + e^{(\alpha (w - W_{refuge}))}}
#'          
#'      Here \eqn{W_{refuge}} marks the body weight at which refuge becomes 
#'      scarcer for prey. \eqn{r} defines the maximum proportion of fish with
#'      access to predation refuge and is always less than or equal to
#'      `max_protect`. \eqn{\alpha} controls the rate at which the 
#'      availability of refuge decreases with increasing body size. It
#'      defaults to a steep slope of 100.
#'
#'      For this method, `method_params` should contain columns named
#'      `prop_protect` and `L_refuge` that give the values for \eqn{r}
#'      and the length at which refuge becomes scarce in cm.
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
#'      `start_L`and `end_L` which contain the starting and ending lengths [cm]
#'      of each size bin and `prop_protect`, the proportion of fish protected
#'      within each corresponding size bin.
#'
#'      3. ``competitive`` - This method is appropriate when refuge density 
#'      data is available for the modelled reef. The refuge density describes 
#'      the distribution of refuges \eqn{(#/m^2)} across defined fish body size
#'      categories. The proportion of fish in size class \eqn{k} with access 
#'      to refuge is given by:
#'
#'      \deqn{R_j(w_p) = \tau\cdot\frac{\eta_k}{\sum_{i}\int_{w_{k-1}}^{w_k}
#'      N_i(w)~dw}  ~~~~~~~ w_p ∈ (~w_{k-1}, w_k~]}
#'      {R_j(w_p) = \tau \eta_k/(\sum_{i}\int_{w_{k-1}}^{w_k}
#'      N_i(w)~dw)  ~~~~~~~ w_p ∈ (~w_{k-1}, w_k~]}
#'
#'      where \eqn{\tau} is the proportion of fish with access to refuge that
#'      are expected to actually utilize  it, \eqn{\eta_k} is the density of
#'      refuges in size range \eqn{(w_{k-1}, w_k]} and
#'       \eqn{\sum_{i}\int_{w_{k-1}}^{w_k} N_i(w)~dw} gives the density
#'       of fish from any group in size range \eqn{(w_{k-1}, w_k]}.
#'      This represents the density of competitors for refuges in
#'      size class \eqn{k}.
#'
#'      For this method, `method_params` should contain columns named
#'      `start_L`and `end_L` which contain the starting and ending lengths [cm]
#'      of each size bin and `refuge_density`, the number of refuges available
#'      in each size bin [no/m^2].
#'
#'  This function checks that the supplied refuge parameters are valid, adds
#'  relevant columns to the `species_params` data frame, and stores refuge
#'  parameters in the `other_params` slot of the `params` object.
#'
#'  Refuge profile parameters can be input in a spreadsheet program and saved
#'  as a .csv file. The data can then be read into R using the command
#'  `read.csv()`.
#'
#' @param params MizerParams object
#' @param method The desired method for setting up benthic refuge,
#'  can be "sigmoidal", "binned", or "competitive"
#' @param method_params A data frame containing values specific to each
#' method for calculating refuge
#' @param refuge_user A vector of logical values indicating whether each
#' functional group uses refuge, TRUE indicates the group uses refuge while
#' false indicates that they do not. Alternatively, these can be included in
#' the `params@species_params` data frame.
#' @param bad_predator A vector of logical values indicating whether hunting is
#' inhibited by refuge for this functional group. FALSE indicates that this
#' species is able to encounter prey within refuge (e.g. eels). Alternatively,
#' these can also be included in the `params@species_params` data frame.
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
#' @concept refuge
#' @family functions for setting parameters
setRefuge <- function(params,
                      method,
                      w_settle = NULL,
                      max_protect = NULL,
                      tau = NULL,
                      method_params,
                      refuge_user = NULL,
                      bad_predator = NULL,...) {

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
    method_options <- c('sigmoidal','binned','competitive')
    if(is.null(method)) {
        stop("You must provide the method to calculate the refuge profile.")
    } else if(!is.element(method, method_options)) {
        stop("Method must be 'sigmoidal','binned', or 'competitive'.")
    }

    # Set default values for parameters used by all methods

    # Minimum size of fish protected by refuges at measured scale
    if(is.null(w_settle)){
        w_settle <- 0.01
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
        tau <- 1
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

    # Check names of method_params for sigmoidal method
    if (refuge_params$method == "sigmoidal") {
        if(!("prop_protect" %in% cnames)) {
            stop("The sigmoidal method parameters dataframe needs a column called
                 'prop_protect' with the proportion of fish protected.")
        } else if(method_params$prop_protect < 0 ||
                  method_params$prop_protect > 1) {
            stop("prop_protect should be a proportion between 0 and 1")
        }
        if(!("L_refuge" %in% cnames)) {
            stop("The sigmoidal method parameters dataframe needs a column called
                 'L_refuge' with the threshhold length (cm) for protected fish.")
        }
        if (is.null(method_params$slope)){ method_params$slope <- 100 }
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
    if (refuge_params$method == "competitive") {
        if(!("start_L" %in% cnames)) {
            stop("The competitive method parameters dataframe needs a 
            column called 'start_L' with the starting lengths (cm) for 
            each size bin.")
        }
        if(!("end_L" %in% cnames)) {
            stop("The competitive method parameters dataframe needs a 
            column called 'end_L' with the end lengths (cm) for
            each size bin.")
        }
        if(!("refuge_density" %in% cnames)) {
            stop("The competitive method parameters dataframe needs a 
            column called 'refuge_density' with the proportion of fish 
            protected for each bin.")
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
            stop("The bad_pred values should be logical.")
        }
        if(length(refuge_user) != no_sp) {
            stop("bad_pred should have a value for every group.")
        }
        params@species_params$bad_pred <- bad_pred
    }

    # Store in params
    params@other_params[['refuge_params']] <- as.data.frame(refuge_params)

    params@other_params[['method_params']] <- as.data.frame(method_params)

    params@species_params$w_settle <- w_settle

    params@time_modified <- lubridate::now()

    return(params)
}

#' Checks unstructured resource parameters and interaction matrix
#'
#' @section Adding unstructured resources:
#'
#'      mizerReef supports two resource spectra that are not size- structured.
#'      Algae are consumed by herbivorous fish, while detritus is consumed by
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
#' @inheritSection getDetritusConsumption Detritus consumption
#' @inheritSection getDetritusProduction Detritus production
#' @inheritSection algal_consumption Algal consumption
#'
#' @param params MizerParams object
#' @param UR_interaction Interaction matrix for unstructured resources
#'  (species x resource)
#'
#' Optional parameters: 
#' @param initial_algae_growth  The initial growth rate of algae in 
#'                              grams/year/m^-2. This value is reset to match
#'                              consumption in the `[reefSteady()]`  function 
#'                              so that steady state abundances match given 
#'                              values.
#'                              
#' @param algae_capacity    The carrying capacity of the system for algal
#'                          biomass in grams per year.
#'                      
#' @param sen_decomp    The proportion of decomposing mass from senescence
#'                      mortality that decomposes to become part of the 
#'                      detritus pool. Defaults to 0.8.
#'                      
#' @param ext_decomp    The proportion of decomposing mass from external
#'                      mortality that decomposes to  become part of the 
#'                      detritus pool. Defaults to 0.2.
#'                      
#' @param initial_d_external    The rate at which detritus biomass sinks from 
#'                              the pelagic zone and becomes part of the 
#'                              detritus pool in grams per year. This value is 
#'                              reset to make up any differences in consumption 
#'                              and production in the `[reefSteady()]` function 
#'                              so that steady state abundances match observed 
#'                              values.
#'
#' @return `setUResourceParams`: MizerParams object with updated unstructured
#'  resource parameters
#' @concept Uresources 
#' @export
setURParams <- function(params,
                        # Preference for resource
                        UR_interaction = NULL, 
                        # Resource Production
                        initial_algae_growth = NULL, 
                        algae_capacity = NULL,
                        sen_decomp = NULL, ext_decomp = NULL, 
                        initial_d_external = NULL) {
    
    # object check ----
        # Check if mizerParams is valid
        assert_that(is(params, "MizerParams"))

        # Find number of species for checks
        no_sp = nrow(params@species_params)
        
    # interaction ----
    # Check if user included in species params
    res_cols <- c('interaction_algae','interaction_detritus')
    if (any(!res_cols %in% names(params@species_params))) {
        # if not, check for provided values
        if(is.null(UR_interaction)){
            stop("You have not provided information on the interaction of
                 functional groups with unstructured resources. You need to 
                 either include and 'interaction_algae' and 
                 'interaction_detritus' columns to species_params or pass
                 a data frame with the same columns to this function.")
        } else {
            # check if matrix has right names
            if(!setequal(names(UR_interaction), res_cols)){
                stop("UR_interaction needs to have columns named 'algae' and
                     'detritus'.")
            }
            # check if detritus and algae interactions are right length
            if(length(UR_interaction$algae) != no_sp) {
                stop("The 'algae' column should have a value for every 
                     functional group.")
            }
            if(length(UR_interaction$detritus) != no_sp){
                stop("The 'detritus' column should have a value for every 
                     functional group.")
            }
            # Check if values are numeric
            if(!all(sapply(UR_interaction, is.numeric))) {
                stop("The entries for algae & detritus interaction should 
                     be numeric.")
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
            
            # Add values as columns to species params data frame
            params@species_params$interaction_algae <- UR_interaction$algae
            params@species_params$interaction_detritus <- UR_interaction$detritus
        }
    }
        
    # set plotting aesthetics ---- 
        params@linecolour["detritus"] <- "plum4"
        params@linecolour["algae"]    <- "olivedrab3"
        params@linetype["detritus"] <- "solid"
        params@linetype["algae"]    <- "solid"
        
    # other parameters ----
        
        ## Production ----
        # Set default algae growth rate
        if(is.null(initial_algae_growth)){ 
            params@other_params$initial_algae_growth <- 2e3
        } else {
            if (!is.numeric(initial_algae_growth)){
                stop("initial_algae_growth should be a numerical value.")
            }
            if (initial_algae_growth <0){
                stop("initial_algae_growth must be non-negative.")
            }
            params@other_params$initial_algae_growth <- initial_algae_growth 
        }
        
        # Set default algae carrying capacity
        if(is.null(algae_capacity)){ 
            params@other_params$algae_capacity <- 2e3
        } else {
            if (!is.numeric(algae_capacity)){
                stop("algae_capacity should be a numerical value.")
            }
            if (algae_capacity < 0){
                stop("algae_capacity must be non-negative.")
            }
            params@other_params$algae_capacity <- algae_capacity 
        }
    
        # Set default proportion of senescence mortality that becomes part of 
        # the detritus pool
        if(is.null(sen_decomp)){ params@other_params$sen_decomp <- 0.8
        } else {
            if (!is.numeric(sen_decomp)){
                stop("sen_decomp should be a numerical value.")
            }
            if (sen_decomp < 0 || sen_decomp > 1){
                stop("sen_decomp must be a proportion between 0 and 1.")
            }
            params@other_params$sen_decomp <- sen_decomp
        }
        
        # Set default proportion of external mortality that becomes part of 
        # the detritus poo
        if(is.null(ext_decomp)){ params@other_params$ext_decomp <- 0.8
        } else {
            if (!is.numeric(ext_decomp)){
                stop("ext_decomp should be a numerical value.")
            }
            if (ext_decomp < 0 || ext_decomp > 1){
                stop("ext_decomp must be a proportion between 0 and 1.")
            }
            params@other_params$ext_decomp <- ext_decomp
        }
    
        # Set default external detritus
        if(is.null(initial_d_external)){ 
            params@other_params$initial_d_external <- 1
        } else { 
            if (!is.numeric(initial_d_external)){
                stop("initial_d_external should be a numerical value.")
            }
            if (initial_d_external < 0){
                stop("initial_d_external must be non-negative.")
            }
            params@other_params$initial_d_external <- initial_d_external 
        }
    
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
#' External mortality is implemented by default in mizerReef. You do not need 
#' to set these parameters. This function should only be used to change default
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
#' @return MizerParams object with updated mortality parameters
#' @concept extmort
#' @export
setExtMortParams <- function(params,
                             ext_mort_params = NULL) {

    # object - Check if mizerParams is valid ----
    assert_that(is(params, "MizerParams"))

    # mort_params - Check if user provided valid mortality parameters ----
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
        # Residual natural mortality
        ext_mort_params$nat_mort  <- 0.2
        # Senescence
        ext_mort_params$sen_prop  <- 0.1
        ext_mort_params$sen_curve <- 0.3
    }

    # Store in params ----
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
#'  mizerReef determines how much fish weight a refuge can hold by converting
#'  user provided fish length bins to weight bins with average length to weight 
#'  conversion parameters \eqn{\bar{a}} and \eqn{\bar{b}}, which describe 
#'  the length to weight relationship for the fish dummies used during data 
#'  collection. The default values are \eqn{\bar{a} = 0.025, \bar{b} = 3}.
#'  Assuming that fish in shelters are neutrally buoyant, the mass of the 
#'  largest fish in each size category is proportional to the volume of 
#'  refuges classified 
#'  in that category. 
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
#'  The mizerReef package provides three methods to define the refuge profile.
#'  
#'  \itemize{
#'
#'  \item **Sigmoidal Method**: \cr
#'  This method is preferred for data-poor reefs or reefs
#'  where the refuge distribution is unknown. It is also ideal for systems
#'  where only one functional group is expected to be utilizing refuge. The
#'  proportion of fish with access to refuge \eqn{ R_j(w_p) } is given by
#'  
#'  \deqn{ R_j(w_p) = \frac{r}
#'                         {1 + e^{ \left( \Delta (w - W_{refuge}) \right)} }}{
#'       R_j(w_p) = r / (1 + e^{( Δ (w - W_{refuge}))})}
#'      
#'  Here \eqn{W_{refuge}} marks the body weight at which refuge becomes 
#'  scarcer for prey. \eqn{r} defines the maximum proportion of fish with
#'  access to predation refuge and is always less than or equal to
#'  `max_protect`. \eqn{\alpha} controls the rate at which the 
#'  availability of refuge decreases with increasing body size. It
#'  defaults to a steep slope of 100.
#'  
#'  For this method, `method_params` should contain columns named
#'  `prop_protect` and `L_refuge` that give the values for \eqn{r}
#'  and the length at which refuge becomes scarce in cm.
#'  
#'  \item **Binned Method**: \cr
#'  This method is appropriate for theoretical applications
#'  and does not rely on empirical data. It sets refuge to a constant
#'  proportion of fish within a given size range. The  proportion of fish
#'  in group \eqn{j} with access to refuge is given by
#'  
#'  \deqn{ R_j(w_p) = r_k ~~~~~~~ w_p ∈ (~w_{k-1}, w_k~] }{
#'       R_j(w_p) = r_k ~~~~~~~ w_p ∈ (~w_{k-1}, w_k~] }
#'       
#'  where \eqn{r_k} is the proportion of fish with access to refuge in
#'  size class \eqn{k}.
#'  
#'  For this method, `method_params` should contain columns named
#'  `start_L`and `end_L` which contain the starting and ending lengths [cm]
#'  of each size bin and `prop_protect`, the proportion of fish protected
#'  within each corresponding size bin.
#'  
#'  \item **Competitive Method**: \cr
#'  This method is appropriate when refuge density 
#'  data is available for the modelled reef. The refuge density describes 
#'  the distribution of refuges \eqn{(no./m^2)} across defined fish body size
#'  categories. The proportion of fish in size class \eqn{k} with access 
#'  to refuge is given by
#'  
#'  \deqn{R_{j}(w_p) = \tau \cdot \frac{ \eta_{k} }
#'                                     { \sum_i \int_{w_{k-1}}^{w-k} N_i(w) \, dw}}{
#'       R_{j}(w_p) = \tau \eta_{k} / 
#'                   ( \sum_i \int_{w_{k-1}}^{w-k} N_i(w) \, dw ) }
#'                            
#'  where \eqn{ \tau } is the proportion of fish with access to refuge that
#'  are expected to actually utilize  it, \eqn{ \eta_{k}} is the density of
#'  refuges in size range \eqn{(w_{k-1}, w_k]} and
#'   \eqn{\sum_{i} \int_{w_{k-1}}^{w_k} N_i(w)~dw} gives the density
#'   of fish from any group in size range \eqn{(w_{k-1}, w_k]}.
#'  This represents the density of competitors for refuges in
#'  size class \eqn{k}.
#'  
#'  For this method, `method_params` should contain columns named
#'  `start_L`and `end_L` which contain the starting and ending lengths [cm]
#'  of each size bin and `refuge_density`, the number of refuges available
#'  in each size bin (no/m^2).
#'      
#'  }
#'
#'  Users can also set a noncomplex reef with no habitat refuge. This option is
#'  convenient for finding steady state parameters. 
#'  
#'  This function checks that the supplied refuge parameters are valid, adds
#'  relevant columns to the `species_params` data frame, and stores refuge
#'  parameters in the `other_params` slot of the `params` object.
#'
#'  Refuge profile parameters can be input in a spreadsheet program and saved
#'  as a .csv file. The data can then be read into R using the command
#'  `read.csv()`.
#'
#' @param params    MizerParams object
#' 
#' @param method    The desired method for setting up benthic refuge, can be 
#'                  "sigmoidal", "binned", "competitive", or "noncomplex"
#'                  
#' @param method_params     A data frame containing values specific to each
#'                          method for calculating refuge
#'                          
#' @param refuge_user   A vector of logical values indicating whether 
#'                      each functional group uses refuge, TRUE indicates the 
#'                      group uses refuge while false indicates that they do 
#'                      not. Must be included here if not in the 
#'                      `params@species_params` data frame.
#'                      
#' @param bad_pred  Optional. A vector of logical values indicating whether 
#'                  hunting is inhibited by refuge for this functional group. 
#'                  FALSE indicates that this species is able to encounter prey 
#'                  within refuge (e.g. eels). Must be included here if not in 
#'                  the `params@species_params` data frame.
#'                  
#' @param satiation Optional. A vector of logical values indicating whether 
#'                  feeding level should be turned on for this functional group.
#'                  Must be included here if not in the 
#'                  `params@species_params` data frame.
#'
#' @param w_settle  Optional. The body weight (g) at which fish settle onto the 
#'                  reef. Fish smaller than this are considered to be larval 
#'                  and thus too small to use predation refuge. Defaults to 0.l 
#'                  grams.
#'                  
#' @param max_protect   Optional. The maximum proportion of fish 
#'                      (any size class) protected by refuge. Defaults to 0.98.
#'                      
#' @param tau   Optional. The proportion of fish with access to refuge that 
#'              actually use it. Defaults to 1.
#'              
#' @param a_bar Optional. The average length to weight conversion value 
#'              \eqn{\bar{a}} and \eqn{\bar{b}} describe the fish dummies used 
#'              when counting refuge holes. \eqn{\bar{a}} defaults to 0.025.
#'              
#' @param b_bar Optional. The average length to weight conversion value 
#'              \eqn{\bar{a}} and \eqn{\bar{b}} describe the fish dummies used 
#'              when counting refuge holes. \eqn{\bar{b}} defaults to 3.
#' @param ... unused
#'
#' @return  A MizerParams object with updated refuge parameters
#' @concept refugeParams
#' @export
setRefuge <- function(params, method, method_params = NULL,
                      # Parameters specific to each group
                      refuge_user = NULL, bad_pred = NULL, 
                      satiation = NULL,
                      # Parameters used by all methods
                      a_bar = NULL, b_bar = NULL,
                      w_settle = NULL, max_protect = NULL, tau = NULL,...) {

    # object check 
        # Check if given mizerParams object is valid
        assert_that(is(params, "MizerParams"))
    
        # Check that a and b parameters are present for all species -
        # needed for l2w conversion
        if (any(!c("a", "b") %in% names(params@species_params))) {
            stop("species_params slot must have columns 'a' and 'b' for ",
                 "length-weight conversion")
        }
        if (anyNA(params@species_params[["a"]]) ||
            anyNA(params@species_params[["b"]])) {
            message("There are NAs in the species_params 
                    columns 'a' and 'b'. You must provide these parameters 
                    for these species if you want them to use refuge.")
        }
    
    # Find number of species for checks
    no_sp = nrow(params@species_params)
    
    # species_params checks
        # Check that refuge_user is logical and the right length
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
        
        # Check that bad_pred is logical and the right length
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
    
        # Check that satiation is logical and the right length
        if(!('satiation' %in% colnames(params@species_params))){
            if(is.null(satiation)){
                stop("You need to provide values for satiation")
            } else if (!is.logical(satiation)) {
                stop("The satiation values should be logical.")
            }
            if(length(satiation) != no_sp) {
                stop("satiation should have a value for every group.")
            }
            params@species_params$satiation <- satiation
        }
    
    # refuge_params set up and checks ----
        # Check if the user provided one of the available refuge profile methods
        method_options <- c('sigmoidal','binned','competitive','noncomplex')
        if(is.null(method)) {
            stop("You must provide the method to calculate the refuge profile.")
        } else if(!is.element(method, method_options)) {
            stop("Method must be 'sigmoidal','binned', 'competitive', 'noncomplex'.")
        }

        # Set default values for parameters used by all methods

        # Minimum weight of fish protected by refuges at measured scale
        if(is.null(w_settle)){
            w_settle <- 0.1
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
            if(!is.numeric(max_protect)) { 
                stop("max_protect should be numeric.") 
            }
            if(max_protect < 0 || max_protect > 1) {
                stop("max_protect should be a proportion between 0 and 1")
            }
        }

        # Proportion of fish with access to refuge that are expected 
        # to utilize it
        if(is.null(tau)){
            tau <- 1
        } else {
            if(!is.numeric(tau)) { stop("tau should be numeric.") }
            if(tau < 0 || tau > 1) {
                stop("tau should be a proportion between 0 and 1")
            }
        }
        
        # a_bar for fish dummies
        if(is.null(a_bar)){
            a_bar <- 0.025
        } else {
            if(!is.numeric(a_bar)) { stop("a_bar should be numeric.") }
            if(a_bar < 0) {
                stop("a_bar must be non-negative.")
            }
        }
        
        # a_bar for fish dummies
        if(is.null(b_bar)){
            b_bar <- 3
        } else {
            if(!is.numeric(a_bar)) { stop("b_bar should be numeric.") }
            if(a_bar < 0) {
                stop("b_bar must be non-negative.")
            }
        }
        
        # Store all values in refuge_params data frame
        refuge_params <- data.frame(method, a_bar, b_bar, 
                                    w_settle, max_protect, tau)
    
    #  method_params set up and checks ----
    if (method != "noncomplex"){
    
        # check if method_params provided
        if (is.null(method_params)) {
            stop("You must provide method specific parameters.")
        }
        
        # check if method_params values are positive and numeric
        if (!is.matrix(method_params)) {
            mp <- as.matrix(method_params)
            # Check that values are numbers
            if (!all(sapply(mp, is.numeric))) {
                stop("The method parameters should be numeric.")
            }
            # Check that all values of refuge method matrix are positive
            if (!all(mp >= 0)) {
                stop("The method parameters must be non-negative.")
            }
        }
        
        # Store column names of method_params for checking
        cnames = colnames(method_params)

        ## Sigmoidal method ----
        if (refuge_params$method == "sigmoidal") {
            # Prop protect
            if(!("prop_protect" %in% cnames)) {
                stop("The sigmoidal method parameters dataframe needs a 
                column called 'prop_protect' with the maximum proportion of 
                fish to be protected.")
            } else if(method_params$prop_protect < 0 ||
                      method_params$prop_protect > 1) {
                        stop("prop_protect should be a proportion between 0 and 1.")
            }
            # L_refuge
            if(!("L_refuge" %in% cnames)) {
                stop("The sigmoidal method parameters dataframe needs a column 
                called 'L_refuge' with the threshhold length (cm) for protected 
                fish.")
            }
            # Slope
            if (is.null(method_params$slope)){ method_params$slope <- 100 } 
        }

        # Binned method ----
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

        # Competitive method ----
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
                column called 'refuge_density' with the density of refuges in 
                each bin in no./square meter.")
            }
            if (!all(method_params$start_L < method_params$end_L)) {
                stop("All bin start lengths must be less than bin end lengths.")
            }
        }
    } else {
        method_params <- as.data.frame(method)
    }
    
    
    # Store in params object --
    params@other_params[['refuge_params']] <- as.data.frame(refuge_params)

    params@other_params[['method_params']] <- as.data.frame(method_params)

    params@time_modified <- lubridate::now()

    return(params)
}

#' Finds the refuge length bins by functional group and stores them
#' params
#'
#' This function is designed to be used after refuge parameters are set by the
#'  [setRefuge()] function. It calculates the proportion of fish 
#' that are in predation refuge for the density-independent sigmoidal and 
#' binned methods. For the competitive method, it finds the indices of fish
#' within the prescribed size bins. These values are used by [reefVulnerable()]
#' to set the vulnerability to predation at each time step.
#' 
#' For all methods, this function calculates the starting and ending body
#' lengths which have access to refuge k. These are calculated with the `a`
#' and `b` parameters specific to each functional group. The lengths are 
#' stored in a data frame called `refuge_lengths` in the `other_params` slot 
#' of the `params` object.
#' 
#' @inheritSection setRefuge Setting the refuge profile
#' 
#' @param params a mizer params object
#' @param ... Unused
#'
#' @return A mizer params object with updated refuge profiles
#' @concept refugeParams
#' @export
getRefuge <- function(params, ...) {
    
    # object - Check if mizerParams is valid ----
    assert_that(is(params, "MizerParams"))
    
    # Extract relevant data from params
    refuge_params <- params@other_params[['refuge_params']]
    method_params <- params@other_params[['method_params']]
    
    # Pull values from params
    w <- params@w
    sp <- params@species_params
    no_w <- length(params@w)
    no_sp <- dim(params@interaction)[1]
    
    # Set parameters used with all methods
    a_bar       <- refuge_params$a_bar
    b_bar       <- refuge_params$b_bar
    w_settle    <- refuge_params$w_settle
    max_protect <- refuge_params$max_protect
    tau         <- refuge_params$tau
    
    # Store which functional groups use refuge
    refuge_user <- sp$refuge_user
    
    # Noncomplex reef ----------------------------------------------------------
    if (refuge_params$method == "noncomplex"){
        
        # Create matrix to store proportions for each species
        refuge <- matrix(0, nrow = no_sp, ncol = no_w)
        rownames(refuge) <- rownames(params@initial_n)
        colnames(refuge) <- colnames(params@initial_n)
        
        # store refuge and bin indices in params object
        params@other_params[['refuge']] <- refuge
        
    # Sigmoidal method ---------------------------------------------------------
    } else if (refuge_params$method == "sigmoidal"){
        
        # Pull slope and proportion of fish to be protected from method_params
        prop_protect <- method_params$prop_protect
        slope <- method_params$slope
        
        # Convert length to weight to determine refuge capacity
        W_refuge <- a_bar * method_params$L_refuge ^ b_bar
        
        # Calculate sigmoid using threshold weights - no organisms smaller 
        # than w_settle or larger than W_refuge can utilize refuge
        denom <- 1 + exp(slope*(w - W_refuge))
        ref <- ifelse(w > w_settle, prop_protect/denom, 0)
        
        # Make sure none of the values are higher than maximum protection allowed
        ref[ref > max_protect] = max_protect
        
        # Create matrix to store proportions for each species
        refuge <- matrix(rep(ref), nrow = no_sp, ncol = no_w, byrow = TRUE)
        rownames(refuge) <- rownames(params@initial_n)
        colnames(refuge) <- colnames(params@initial_n)
        
        # Account for species that don't utilize refuge
        refuge <- refuge_user*refuge
        
        # store refuge in params object
        params@other_params[['refuge']] <- refuge
        
        # Find L_refuge by species & store in data frame
        L_refuge.i <- (W_refuge / sp[["a"]])^(1 / sp[["b"]])
        refuge_lengths <- data.frame(sp$species, L_refuge.i)
        params@other_params[['refuge_lengths']] <- refuge_lengths
        
        # Save time parameters were modified
        params@time_modified <- lubridate::now()
        
    # Binned method ------------------------------------------------------------
    } else if (refuge_params$method == "binned") {
        
        # Initialize storage
        ref       <- rep(0, no_w)
        start_l.i <- list(1)
        end_l.i   <- list(1)
        bin.id   <- list (1)
        no_bins   <- nrow(method_params)
        
        # Loop through each refuge bin
        for (k in 1:no_bins) {
            
            # Calculate start and end of weight bins for a dummy fish
            start_w <- a_bar * method_params$start_L[[k]] ^ b_bar
            end_w <- a_bar * method_params$end_L[[k]] ^ b_bar
            
            # Set threshold weight - no organisms smaller than w_settle
            start_w[start_w < w_settle] <- w_settle
            
            # Gives indices of fish in size range to protect
            bin.id[[k]] <- which(w >= start_w & w <= end_w)
            
            # Refuge
            ref[bin.id[[k]]] = method_params$prop_protect[k]
            
            # Calculate length bins for each species
            start_l.i[[k]] <- (start_w / sp[["a"]])^(1 / sp[["b"]])
            names(start_l.i)[[k]] <- c(paste("start",k,sep = ""))
            end_l.i[[k]]   <- (end_w   / sp[["a"]])^(1 / sp[["b"]])
            names(end_l.i)[[k]] <- c(paste("end",k,sep = ""))
        }
        
        # Create matrix to store proportions for each species
        refuge <- matrix(rep(ref), nrow = no_sp, ncol = no_w, byrow = TRUE)
        rownames(refuge) <- rownames(params@initial_n)
        colnames(refuge) <- colnames(params@initial_n)
        
        # Make sure none of the values are higher than maximum protection allowed
        refuge[refuge > max_protect] = max_protect
        
        # Account for species that don't utilize refuge
        refuge <- refuge_user*refuge
        
        # store refuge and bin indices in params object
        params@other_params[['refuge']] <- refuge
        params@other_params[['bin.id']] <- bin.id
        
        # store length bins by functional group in params object
        start_l.i <- t(do.call(rbind, start_l.i))
        end_l.i  <- t(do.call(rbind, end_l.i))
        refuge_lengths <- cbind(start_l.i, end_l.i)
        row.names(refuge_lengths) <- sp$species
        params@other_params[['refuge_lengths']] <- refuge_lengths
        
        # Save time parameters were modified
        params@time_modified <- lubridate::now()
        
    ## Competitive method ------------------------------------------------------
    } else if (refuge_params$method == "competitive") {
        
        # Initialize empty list to hold number of competitors for each bin
        competitor_density = numeric(nrow(method_params))
        
        # Empty list to hold indices of fish protected by each bin
        bin.id = list(1)
        start_l.i <- list(1)
        end_l.i   <- list(1)
        
        # Loop through each refuge bin
        for (k in 1:nrow(method_params)) {
            
            # Calculate start and end of weight bin k
            start_w <- a_bar * method_params$start_L[[k]] ^ b_bar
            end_w <- a_bar * method_params$end_L[[k]] ^ b_bar
            
            # No organisms smaller than w_settle can use refuge
            start_w[start_w < w_settle] <- w_settle
            
            # Find indices of fish within size bin k
            bin.id[[k]] <- which(params@w >= start_w & params@w <= end_w)
            
            # Calculate length bins for each species
            start_l.i[[k]] <- (start_w / sp[["a"]])^(1 / sp[["b"]])
            names(start_l.i)[[k]] <- c(paste("start",k,sep = ""))
            end_l.i[[k]]   <- (end_w   / sp[["a"]])^(1 / sp[["b"]])
            names(end_l.i)[[k]] <- c(paste("end",k,sep = ""))
        }
        
        # Store indices of each bin
        params@other_params[['bin.id']] <- bin.id
        
        # Store length bins by functional group in a data frame
        start_l.i <- t(do.call(rbind, start_l.i))
        end_l.i  <- t(do.call(rbind, end_l.i))
        refuge_lengths <- cbind(start_l.i, end_l.i)
        row.names(refuge_lengths) <- sp$species
        params@other_params[['refuge_lengths']] <- refuge_lengths
        
        # Save time parameters were modified
        params@time_modified <- lubridate::now()
    }
    return(params)
}

#' Change the refuge parameters for a model
#'
#' This is a wrapper function for the [setRefuge()] and [getRefuge()] 
#' functions that allows users to easily change refuge parameters on an
#' existing mizer model. This will take it out of steady state, and users
#' should run [reefSteady()] after this function to return to steady state. 
#' 
#' @inheritSection setRefuge Setting the refuge profile
#' 
#' @param params a mizer object
#' 
#' @param new_method    The new method to be used for setting the refuge 
#'                      profile. Options are "sigmoidal", "binned", 
#'                      "competitive", or "noncomplex". If no method is
#'                      provided, this defaults to the same method as is
#'                      currently being used in the simulation. 
#'                  
#' @param new_method_params  A data frame containing values specific to each
#'                          method for calculating refuge. Only necessary if 
#'                          changing methods. 
#'                          
#' @param new_L_refuge  To be used with "sigmoidal" method only. The new value
#'                      for the length at which refuge becomes scarce in cm.
#'                      
#' @param new_prop_protect  To be used with "sigmoidal" method only. The new 
#'                          value for the maximum proportion of fish to protect.
#'                          
#' @param scale_bin To be use with the "binned" method only. A number or vector 
#'                  of numbers to multiply the `prop_protect` values by for 
#'                  each size bin. Changes the proportion of fish protected.
#'                          
#' @param ... Unused
#'
#' @return A mizer object with updated refuge profiles
#' @concept refugeParams
#' @export
newRefuge <- function(params,
                      # Fully changing method
                      new_method = NULL, new_method_params = NULL,
                      # Sigmoidal - changing refuge length prop protect
                      new_L_refuge = NULL, new_prop_protect = NULL,
                      # Binned - scaling bin prop 
                      scale_bin = NULL, ...) {
    
    # Check that the user provided at least one new input
        inputs <- list(new_method, new_method_params, 
                    new_L_refuge, new_prop_protect,
                    scale_bin)
        if (all(sapply(inputs, is.null))) {
            stop("Error: At least one input must be provided.")
        }
    
      # object check ----  
        # Check if given mizerParams object is valid
        assert_that(is(params, "MizerParams"))
        
        # Extract relevant data from params
        refuge_params <- params@other_params[['refuge_params']]
        method_params <- params@other_params[['method_params']]
    
    # method checks ----
        # Check if the user provided one of the available refuge profile methods
        m_opt <- c('sigmoidal','binned','competitive','noncomplex')
        if(is.null(new_method)) {
            warning("Since you did not provide a new method for determinging 
                    the refuge profile, I will use the one currently stored in 
                    params@other_params$refuge_params$method.")
            # If user did not provide a method, use old one
            new_method <- refuge_params$method
        # If user did provide a method, check that it's one of the options
        } else if(!is.element(new_method, m_opt)) {
            stop("Method must be 'sigmoidal','binned', 'competitive', 'noncomplex'.")
        }
    
    #  method_params checks ----
        if (new_method != "noncomplex"){
            # check if method_params provided
            if (is.null(new_method_params)) {
                # Return an error if they are trying to switch to methods
                if(refuge_params$method != new_method){
                    stop("You must provide a new method_params data frame to
                         switch methods.") 
                }
                # Check method ----
                ## Sigmoidal ----
                if (new_method == "sigmoidal") {
                    # Check if new L or new_prop given
                    if (all(is.null(new_L_refuge), is.null(new_prop_protect))){
                       stop("You must provide either a new L_refuge, a new
                             prop_protect, or a new method_params data frame.") 
                    }
                    # If new L_refuge given, check that it's positive
                    if(!is.null(new_L_refuge)){
                        if(new_L_refuge < 0) {
                            stop("new_L_refuge must be non-negative.")
                        }
                    } else { new_L_refuge <- method_params$L_refuge }
                    # If new prop_protect given, check that it's between 0 and 1
                    if(!is.null(new_prop_protect)){
                        if(method_params$prop_protect < 0 ||
                           method_params$prop_protect > 1) {
                        stop("prop_protect should be a proportion between 0 and 1.")
                        }
                    } else { new_prop_protect <- method_params$prop_protect }
                    
                    new_mp <- method_params
                    new_mp$L_refuge <- new_L_refuge
                    new_mp$prop_protect = new_prop_protect
                    
                    ## Binned or Competitive ----
                } else if (new_method == "binned" || new_method == "competitive") {
                    # Find number of bins used in old method
                    no_bins <- length(method_params)
                    if(is.null(scale_bin)){
                        stop("You must provide either a value or vector of values
                             for scale_bin or a new method_params data frame.") 
                    }
                    # Make sure scaling vector is the right length
                    if(!(length(scale_bin) == 1) || 
                       !(length(scale_bin) == no_bins)){
                        stop("scale_bin must have length 1 or have an value
                             for every bin.")
                    }
                    # Check that scale_bin_prop is positive
                    if(scale_bin < 0) {
                        stop("scale_bin must be non-negative.")
                    }
                    # Calculate new bins
                    if (new_method == "binned"){ 
                        new_mp <- method_params
                        new_mp$prop_protect <- scale_bin * new_mp$prop_protect
                        
                    } else {
                        new_mp <- method_params
                        new_mp$refuge_density <- scale_bin * new_mp$refuge_density
                    }
                }
            } else { new_mp <- new_method_params }
        }
        
    #### Update parameters ----
        params <- setRefuge(params = params, 
                            method = new_method,
                            method_params = new_mp)
        
        params <- getRefuge(params = params)
}





#' Set algae parameters for mizerReef
#'
#' @section Algae as an unstructured resource:
#'
#'      mizerReef supports algae as a non-size-structured resource,
#'      consumed primarily by herbivorous fish. This function sets
#'      the initial growth rate, system carrying capacity, and
#'      interaction strengths for algae, allowing for flexible
#'      diet preferences.
#'
#'      The interaction strength (\eqn{\theta_{i,algae}}) for each
#'      species \eqn{i} determines how strongly that group feeds on
#'      algae. This can be set via the `interaction_algae` column in
#'      the species parameter data frame, or directly via the
#'      `UR_interaction` argument. If neither is provided, all
#'      interaction strengths are set to zero and a warning is issued.
#'
#'      The initial growth rate (`algae_growth_initial`) and carrying
#'      capacity (`algae_capacity`) control the baseline production
#'      and maximum standing stock of algae, respectively. These values
#'      may be reset by [reefSteady()] to ensure steady state abundances
#'      match observed or target values.
#'
#'      Carrying capacity can be toggled with `use_UR_cc`. When enabled,
#'      algae biomass will be limited by the specified capacity.
#'
#'      Note: Interaction with size-structured resources, such as plankton,
#'      is set with the resource_interaction column of the species parameters
#'      dataframe.
#'
#' @inheritSection algae_consumption Algae consumption
#'
#' @param params A `MizerParams` object.
#' @param algae_growth_initial Numeric. The initial growth rate of algae
#'                             in grams/m^2/year. Default is 2e3.
#' @param algae_capacity Numeric. Carrying capacity for algae biomass in
#'                       grams per year. Default is 1e4.
#' @param UR_interaction Optional. A named list or array with one or more
#'                       resource interaction vectors (e.g. interaction_algae,
#'                       interaction_detritus, interaction_sponge), each of
#'                       length equal to the number of species. If NULL, will
#'                       use columns in species_params or set to zero. All
#'                       values must be numeric and between 0 and 1.
#' @param use_UR_cc Logical. Whether to implement a carrying capacity for
#'                  all unstructured resources. Default is FALSE. This flag
#'                  is stored in the other_params slot.
#' @param algae_colour Character. Colour to use for algae in plots.
#'                     Default is "darkseagreen3".
#' @details All algae-related parameters (growth rate, capacity) are stored in
#'          the `algae_params` slot of the params object. Resource interaction
#'          strengths are set in the species_params data frame. This function
#'          supports flexible multi-resource interaction via the
#'          UR_interaction argument.
#' @return A `MizerParams` object with updated algae parameters (in `algae_params` slot).
#' @concept algae
#' @export
setAlgaeParams <- function(params,
                           algae_growth_initial = NULL,
                           algae_capacity = NULL,
                           UR_interaction = NULL,
                           use_UR_cc = FALSE,
                           algae_colour = "darkseagreen3") {
    assert_that(is(params, "MizerParams"))
    no_sp <- nrow(params@species_params)
    assert_that(is.flag(use_UR_cc))
    params@other_params$use_UR_cc <- use_UR_cc

    # Ensure algae_params slot exists
    if (is.null(params@algae_params)) params@algae_params <- list()

    # Interaction setup (array for all URs)
    if (is.null(UR_interaction)) {
        if (!("interaction_algae" %in% names(params@species_params))) {
            params@species_params$interaction_algae <- rep(0, no_sp)
            warning("You have not provided any values for interaction_algae, so no species feed on algae.")
        }
    } else {
        # UR_interaction should be a named list/array with possible multiple resources
        if (!is.list(UR_interaction) && !is.array(UR_interaction)) {
            stop("UR_interaction must be a named list or array with resource columns
                 (e.g. interaction_algae, interaction_detritus, interaction_sponge)")
        }
        for (ur in names(UR_interaction)) {
            vals <- UR_interaction[[ur]]
            if (length(vals) != no_sp) stop(paste0(ur, " must have a value for every species."))
            if (!is.numeric(vals)) stop(paste0(ur, " must be numeric."))
            if (any(vals < 0) || any(vals > 1)) stop(paste0(ur, " values must be between 0 and 1."))
            params@species_params[[ur]] <- vals
        }
    }

    # Store algae params in slot with checks
    # algae_growth_initial check
    if (is.null(algae_growth_initial)) {
        params@algae_params$algae_growth <- 2e3
    } else {
        if (!is.numeric(algae_growth_initial)) stop("algae_growth_initial should be a numerical value.")
        if (algae_growth_initial < 0) stop("algae_growth_initial must be non-negative.")
        params@algae_params$algae_growth <- algae_growth_initial
    }

    # algae_capacity check
    if (is.null(algae_capacity)) {
        params@algae_params$algae_capacity <- 1
    } else {
        if (!is.numeric(algae_capacity)) stop("algae_capacity should be a numerical value.")
        if (algae_capacity < 0) stop("algae_capacity must be non-negative.")
        params@algae_params$algae_capacity <- algae_capacity
    }

    # Colour check (for algae_colour)
    if (!is.null(algae_colour)) {
        if (!is.character(algae_colour) || length(algae_colour) != 1 || is.na(grDevices::col2rgb(algae_colour)[1])) {
            stop("algae_colour must be a valid color name or hex code.")
        }
    }
    params@linecolour["algae"] <- algae_colour
    params@linetype["algae"] <- "solid"
    params@time_modified <- lubridate::now()
    return(params)
}

#' Set detritus parameters for mizerReef
#'
#' @section Detritus as an unstructured resource:
#'
#'      Detritus in mizerReef is modeled as a non-size-structured resource
#'      produced by the decomposition of organic materials. Detritus is consumed
#'      by detritivores and benthic invertebrates. This function sets the
#'      carrying capacity, decomposition proportions, and interaction strengths
#'      for detritus, supporting flexible diet preferences.
#'
#'      The interaction strength (\eqn{\theta_{i,detritus}}) for each species \eqn{i}
#'      determines how strongly that group feeds on detritus. This can be set via the
#'      `interaction_detritus` column in the species parameter data frame, or directly
#'      via the `UR_interaction` argument. If neither is provided, all interaction
#'      strengths are set to zero and a warning is issued.
#'
#'      The carrying capacity (`detritus_capacity`) limits the maximum standing
#'      stock of detritus.
#'
#'      The proportions of decomposing mass from senescence (`sen_decomp`) and
#'      external mortality (`ext_decomp`) that become detritus are set here,
#'      with typical defaults of 0.8 and 0.2, respectively. The rate at which
#'      detritus sinks from the pelagic zone (`d_external_initial`) is also set,
#'      controlling detritus input from external sources. These values may be
#'      reset by [reefSteady()] to ensure steady state abundances match observed
#'      or target values.
#'
#'      Carrying capacity can be toggled with `use_UR_cc`. When enabled, detritus
#'      biomass will be limited by the specified capacity.
#'
#'      Note: Interaction with size-structured resources, such as plankton, is
#'      set with the resource_interaction column of the species parameters dataframe.
#'
#' @inheritSection getDetritusProduction Detritus production
#' @inheritSection getDetritusConsumption Detritus consumption
#'
#' @param params A `MizerParams` object.
#' @param detritus_capacity Numeric. Carrying capacity for detritus biomass in grams
#'                          per year. Default is 1e4.
#' @param sen_decomp Numeric. Proportion of decomposing mass from senescence mortality
#'                   that becomes detritus. Default is 0.8.
#' @param ext_decomp Numeric. Proportion of decomposing mass from external mortality
#'                   that becomes detritus. Default is 0.2.
#' @param external Numeric. Rate at which detritus biomass sinks from the pelagic
#'                 zone (grams per year). Default is 1.
#' @param UR_interaction Optional. A named list or array with one or more resource interaction
#'                       vectors (e.g. interaction_algae, interaction_detritus, interaction_sponge),
#'                       each of length equal to the number of species. If NULL, will use columns
#'                       in species_params or set to zero. All values must be numeric and between 0 and 1.
#' @param use_UR_cc Logical. Whether to implement a carrying capacity for detritus. Default
#'                  is FALSE. This flag is stored in the other_params slot.
#' @param detritus_colour Character. Colour to use for detritus in plots. Default is "plum4".
#'
#' @details All detritus-related parameters (capacity, decomposition rates, external input) are
#'          stored in the `detritus_params` slot of the params object. Resource interaction
#'          strengths are set in the species_params data frame. This function supports flexible
#'          multi-resource interaction via the UR_interaction argument.
#' @return A `MizerParams` object with updated detritus parameters (in `detritus_params` slot).
#' @concept detritus
#' @export
setDetritusParams <- function(params,
                              detritus_capacity = NULL,
                              sen_decomp = NULL,
                              ext_decomp = NULL,
                              external = NULL,
                              UR_interaction = NULL,
                              use_UR_cc = FALSE,
                              detritus_colour = "plum4") {
    assert_that(is(params, "MizerParams"))
    no_sp <- nrow(params@species_params)
    assert_that(is.flag(use_UR_cc))
    params@other_params$use_UR_cc <- use_UR_cc

    # Ensure detritus_params slot exists
    if (is.null(params@detritus_params)) params@detritus_params <- list()

    # Interaction setup (array for all URs)
    if (is.null(UR_interaction)) {
        if (!("interaction_detritus" %in% names(params@species_params))) {
            params@species_params$interaction_detritus <- rep(0, no_sp)
            warning("You have not provided any values for interaction_detritus, so no species currently feed on detritus.")
        }
        ur_names <- grep("^interaction_", names(params@species_params), value = TRUE)
        for (ur in ur_names) {
            if (is.null(params@species_params[[ur]])) {
                params@species_params[[ur]] <- rep(0, no_sp)
            }
        }
    } else {
        if (!is.list(UR_interaction) && !is.array(UR_interaction)) {
            stop("UR_interaction must be a named list or array with resource columns (e.g. interaction_algae, interaction_detritus, interaction_sponge)")
        }
        for (ur in names(UR_interaction)) {
            vals <- UR_interaction[[ur]]
            if (length(vals) != no_sp) stop(paste0(ur, " must have a value for every species."))
            if (!is.numeric(vals)) stop(paste0(ur, " must be numeric."))
            if (any(vals < 0) || any(vals > 1)) stop(paste0(ur, " values must be between 0 and 1."))
            params@species_params[[ur]] <- vals
        }
    }

    # Store detritus params in slot with checks
    params@detritus_params$detritus_capacity <- ifelse(is.null(detritus_capacity), 1, detritus_capacity)

    # sen_decomp check
    if (is.null(sen_decomp)) {
        params@detritus_params$sen_decomp <- 0.8
    } else {
        if (!is.numeric(sen_decomp)) stop("sen_decomp should be a numerical value.")
        if (sen_decomp < 0 || sen_decomp > 1) stop("sen_decomp must be a proportion between 0 and 1.")
        params@detritus_params$sen_decomp <- sen_decomp
    }

    # ext_decomp check
    if (is.null(ext_decomp)) {
        params@detritus_params$ext_decomp <- 0.2
    } else {
        if (!is.numeric(ext_decomp)) stop("ext_decomp should be a numerical value.")
        if (ext_decomp < 0 || ext_decomp > 1) stop("ext_decomp must be a proportion between 0 and 1.")
        params@detritus_params$ext_decomp <- ext_decomp
    }

    params@detritus_params$external <- ifelse(is.null(external), 1, external)

    # Colour check (for detritus_colour)
    if (!is.null(detritus_colour)) {
        if (!is.character(detritus_colour) || length(detritus_colour) != 1 || is.na(grDevices::col2rgb(detritus_colour)[1])) {
            stop("detritus_colour must be a valid color name or hex code.")
        }
    }
    params@linecolour["detritus"] <- detritus_colour
    params@linetype["detritus"] <- "solid"
    params@time_modified <- lubridate::now()
    return(params)
}

#' Set the parameters for external mortality
#'
#' @section External mortality in mizerReef:
#'
#'      External mortality in mizerReef includes residual natural mortality
#'      and senescence mortality, representing background sources of death
#'      not explicitly modeled (e.g., predation, disease, aging). This function
#'      allows you to set or override the default rates and scaling for these
#'      processes.
#'
#'      Residual natural mortality is modeled as an allometric function of body
#'      size, decreasing with increasing size. Senescence mortality is modeled
#'      as a power function of relative size, allowing flexible control over
#'      the rate and scaling of age-related death.
#'
#'      If no parameters are provided, defaults are used: nat_mort = 0.2,
#'      sen_prop = 0.1, sen_curve = 0.3. All parameters must be numeric
#'      and non-negative. The function checks for required columns and
#'      valid values if a custom parameter set is provided.
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
#'                              \frac{\log_{10}(w)}{\log_{10}(w_{max.i})}
#'                              \right)^{p_{sen}}}
#'           {\mu_{sen.i}(w) = k_{sen}
#'           (\log_{10}(w)/\log_{10}(w_{max.i}))^{p_{sen}}}
#'
#'      where \eqn{k_{sen}} is the rate of senescence mortality, \eqn{p_{sen}}
#'      defines the slope of the senescence curve, \eqn{w_{max.i}} maximum body
#'      size of group \eqn{i} in grams.
#'
#' @param params A `MizerParams` object.
#' @param ext_mort_params Optional. A named list or matrix with columns/names
#'                        'nat_mort', 'sen_prop', and 'sen_curve'. If NULL,
#'                        defaults are used. All values must be numeric and
#'                        non-negative.
#' @return A `MizerParams` object with updated external mortality parameters
#'         (stored in `other_params$ext_mort_params`).
#' @concept extmort
#' @export
setExtMortParams <- function(params,
                             ext_mort_params = NULL) {
    # object - Check if mizerParams is valid
    assert_that(is(params, "MizerParams"))

    # mort_params - Check if user provided valid mortality parameters
    if (!is.null(ext_mort_params)) {
        if (!all(sapply(ext_mort_params, is.numeric))) {
            stop("The external mortality parameters should be numeric.")
        }
        if (!is.matrix(ext_mort_params)) {
            ext_mort_params <- as.matrix(ext_mort_params)
            if (!all(ext_mort_params >= 0)) {
                stop("The external mortality parameters should be
                     nonnegative")
            }
        }
        if (!("nat_mort" %in% colnames(ext_mort_params))) {
            stop("The external mortality parameters dataframe needs a column
                 called 'nat_mort' with the residual natural mortality rate.")
        }
        if (!("sen_prop" %in% colnames(ext_mort_params))) {
            stop("The external mortality parameters dataframe needs a column
                 called 'sen_prop' with the external mortality rate.")
        }
        if (!("sen_curve" %in% colnames(ext_mort_params))) {
            stop("The external mortality parameters dataframe needs a column
                 called 'sen_curve' with the exponent for the external
                 curve.")
        }
    } else {
        ext_mort_params <- vector("list", 3)
        names(ext_mort_params) <- c("nat_mort", "sen_prop", "sen_curve")
        # Residual natural mortality
        ext_mort_params$nat_mort <- 0.2
        # Senescence
        ext_mort_params$sen_prop <- 0.1
        ext_mort_params$sen_curve <- 0.3
    }

    # Store in params
    params@other_params[["ext_mort_params"]] <- ext_mort_params

    params@time_modified <- lubridate::now()

    return(params)
}

#' Set the refuge profile parameters
#'
#' Sets and validates refuge parameters in the mizerParams object, allowing
#' flexible calibration to match ecological data and measurement approaches.
#' Supports multiple methods for defining refuge profiles and options for
#' calibrating bin boundaries using either species-specific or dummy fish
#' length-weight parameters.
#'
#' Refuge profiles account for the protective behavior of prey living in
#' high-complexity environments (e.g. coral reefs) with access to predation
#' refuge. The refuge profile defines the proportion of fish within user-defined
#' length bins that are protected from being encountered by a predator.
#'
#' A unique refuge profile is generated for each predator group × prey group ×
#' prey size combination based on the given refuge profile parameters and four
#' values from `params@species_params`: length-weight conversion values `a` and
#' `b`, `refuge_user` (TRUE for groups that utilize predation refuge), and
#' `blocked_pred` (FALSE for predator groups whose body shape or predatory
#' strategy allow them to access fish within refuge, e.g. eels).
#'
#' The maximum proportion of fish protected by refuge in any size class is set
#' by `max_protect` to ensure some food is always available to predators.
#'
#' The refuge profile is used when calculating the food encounter rate in
#' [reefEncounter()] and the predation mortality rate in [reefPredMort()]. Its
#' entries are dimensionless values between 0 and 1, representing the proportion
#' of fish in the corresponding prey and size categories that are hidden within
#' refuge and thus cannot be encountered by predators. If no refuge is available,
#' predator-prey interactions are determined entirely by size-preference.
#'
#' @section Defining refuge threshholds:
#'
#'      Defining the refuge profile for a given system is nontrivial. There are many
#'      ways to measure and define predation refuge on reefs. The way you parameterise
#'      the model should reflect your data collection method and ecological context.
#'
#'      The `use_dummy_fish_bins` argument determines how refuge bin boundaries are
#'      set:
#'
#'      - If TRUE (default), bin boundaries are determined by weight, calculated
#'        using dummy fish length-weight parameters (`a_bar`, `b_bar`). All species
#'        share the same weight boundaries for bins, but the lengths of fish in each
#'        bin are calculated for each species using their own length-weight
#'        parameters (`a`, `b`). This means fish in the same bin may be of very
#'        different lengths depending on species. This is appropriate when refuge
#'        capacity is measured using dummy fish of known size and is the
#'        preferred method.
#'
#'      - If FALSE, bin boundaries are determined by length, so all species in a bin
#'        have the same length, but their weights are calculated using their own
#'        length-weight parameters (`a`, `b`). This means fish in the same bin may
#'        have very different weights depending on species. This is appropriate when
#'        refuge hole entrances are measured for each species.
#'
#'  This choice applies to all refuge methods (sigmoidal, binned, competitive)
#'  and is stored in `params@refuge_params$use_dummy_fish_bins`. Select the
#'  option that matches your measurement approach and ecological realism.
#'
#' @section Setting the refuge profile:
#'
#'      The mizerReef package provides three methods to define the refuge profile.
#'
#'      \itemize{
#'
#'      \item **Sigmoidal Method**: \cr
#'
#'          This method is preferred for data-poor reefs or reefs where the refuge
#'          distribution is unknown. It is also ideal for systems where only one
#'          species is expected to be utilizing refuge. The sigmoidal method defines
#'          a smooth transition in refuge availability around a threshold body size.
#'
#'          The threshold for refuge can be set in two ways, depending on how your
#'          refuge data was collected:
#'
#'          - If `use_dummy_fish_bins = FALSE`, the threshold weight is calculated as:
#'                  \deqn{ W_{i.refuge} = a_i \cdot L_{refuge}^{b_i} }
#'           where \eqn{a_i} and \eqn{b_i} are the length-weight parameters for species i.
#'
#'          - If `use_dummy_fish_bins = TRUE`, the weight threshold is:
#'
#'                  \deqn{ W_{refuge} = a_{bar} \cdot L_{refuge}^{b_{bar}} }
#'
#'           The proportion of fish with access to refuge is then given by:
#'                  \deqn{ R_{j}(w_p) = \frac{r}{1 + e^{\Delta(w - W_{refuge})}} }
#'          where $r$ is the maximum proportion protected, $\Delta$ is the slope,
#'          $w$ is body weight, and $W_{refuge}$ is the threshold (species-specific or dummy fish)
#'
#'          For this method, `method_params` should contain columns named
#'          `prop_protect` and `L_refuge` that give the values for \eqn{r}
#'          and the length at which refuge becomes scarce in cm.
#'
#'      \item **Binned Method**: \cr
#'
#'          This method is appropriate for theoretical applications
#'          and does not rely on empirical data. It sets refuge to a constant
#'          proportion of fish within a given size range. The proportion of fish
#'          in group \eqn{j} with access to refuge is given by
#'
#'          \deqn{ R_j(w_p) = r_k ~~~~~~~ w_p ∈ (~w_{k-1}, w_k~] }{
#'               R_j(w_p) = r_k ~~~~~~~ w_p ∈ (~w_{k-1}, w_k~] }
#'
#'          where \eqn{r_k} is the proportion of fish with access to refuge in
#'          size class \eqn{k}.
#'
#'          For this method, `method_params` should contain columns named
#'          `start_L` and `end_L` which contain the starting and ending lengths [cm]
#'          of each size bin and `prop_protect`, the proportion of fish protected
#'          within each corresponding size bin.
#'
#'      \item **Competitive Method**: \cr
#'          This method is appropriate when refuge density data is available for
#'          the modelled reef. The refuge density describes the distribution of
#'          refuges \eqn{(no./m^2)} across predefined fish body size categories.
#'          The proportion of fish in size class \eqn{k} with access to refuge
#'          is given by
#'
#'          \deqn{R_{j}(w_p) = \tau \cdot \frac{ \eta_{k} }
#'                                             { \sum_i \int_{w_{k-1}}^{w-k} N_i(w) \, dw}}{
#'               R_{j}(w_p) = \tau \eta_{k} /
#'                           ( \sum_i \int_{w_{k-1}}^{w-k} N_i(w) \, dw ) }
#'
#'          where \eqn{ \tau } is the proportion of fish with access to refuge that
#'          are expected to actually utilize  it, \eqn{ \eta_{k}} is the density of
#'          refuges in size range \eqn{(w_{k-1}, w_k]} and
#'          \eqn{\sum_{i} \int_{w_{k-1}}^{w_k} N_i(w)~dw} gives the density
#'          of fish from any group in size range \eqn{(w_{k-1}, w_k]}.
#'          This represents the density of competitors for refuges in
#'          size class \eqn{k}.
#'
#'          For this method, `method_params` should contain columns named
#'          `start_L`and `end_L` which contain the starting and ending lengths [cm]
#'          of each size bin and `refuge_density`, the number of refuges available
#'          in each size bin (no/m^2).
#'
#'      }
#'
#'  Users can also set a noncomplex reef with no habitat refuge. This option is
#'  convenient for finding steady state parameters and is the default when
#'  no parameters are provided.
#'
#'  This function checks that the supplied refuge parameters are valid, adds
#'  relevant columns to the `species_params` data frame, and stores refuge
#'  parameters in the `refuge_params` slot of the `params` object.
#'
#'  Refuge profile parameters can be input in a spreadsheet program and saved
#'  as a .csv file. The data can then be read into R using the command
#'  `read.csv()`.
#'
#' @param params    MizerParams object
#' @param method    Character. The method for setting up benthic refuge.
#'                 One of "sigmoidal", "binned", "competitive", or
#'                 "noncomplex". See Details. **Required.**
#' @param method_params    Data frame or named list. Specifies parameters required for the chosen method:
#'   - For "sigmoidal": must include `L_refuge` (numeric, length at which refuge becomes scarce, cm; **no default**) and `prop_protect` (numeric, max proportion protected, default: 0.98).
#'   - For "binned": must include `start_L` (numeric, start length, cm; **no default**), `end_L` (numeric, end length, cm; **no default**), and `prop_protect` (numeric, proportion protected, default: 0.98).
#'   - For "competitive": must include `start_L`, `end_L` (as above), and `refuge_density` (numeric, refuges per size bin, no/m^2; **no default**).
#'   - For "noncomplex": no parameters required.
#'
## Choice of refuge bin calibration
#'
#' @param use_dummy_fish_bins Logical. Controls how refuge bin boundaries and thresholds are calculated for sigmoidal and binned methods:
#'   - TRUE: Use each species' own length-weight parameters (`a`, `b`).
#'   - FALSE (default): Use dummy fish parameters (`a_bar`, `b_bar`).
#'   Set according to your data collection method. The setting is stored in
#'   `params@refuge_params$use_dummy_fish_bins` and used by [getRefuge()].
#'
#' @param refuge_user   Logical vector (length = number of species). Indicates which groups use refuge.
#'                      If not present in `species_params`, must be provided. Defaults to FALSE.
#'
#' @param blocked_pred Optional. Logical vector (length = number of species). Indicates whether the predator is blocked by refuge for this species. TRUE means hunting is blocked by refuge; FALSE means the species can access prey within refuge (e.g. eels). Defaults to FALSE.
#'
#' @param satiation Logical vector (length = number of species). Indicates which groups are subject to
#'                  satiation. If not provided, defaults are set automatically: FALSE for carnivores
#'                  (species that eat other species, i.e. row sum of interaction matrix > 0), TRUE
#'                  for pure resource consumers (species that do not eat other species but have
#'                  positive resource interaction). A warning is issued if defaults are used.
#'
#' @param a_bar Numeric. Length-weight conversion parameter for dummy fish. **Default:** 0.025.
#'              If any species is missing an 'a' parameter, the value of a_bar is used for that
#'              species and a warning is issued.
#'
#' @param b_bar Numeric. Length-weight exponent for dummy fish. **Default:** 3.
#'              If any species is missing a 'b' parameter, the value of b_bar
#'              is used for that species and a warning is issued.
#'
#' @param w_settle Numeric. Minimum weight of fish protected by refuges at measured scale (grams). **Default:** 0.1.
#'
#' @param max_protect Numeric. Maximum proportion of fish protected by refuge (0–1). **Default:** 0.98.
#'
#' @param tau Numeric. Proportion of fish with access to refuge expected to utilize it (0–1). **Default:** 1.
#'
#' @param ... Unused.
#'
#' @seealso
#'   [reefEncounter()], [reefPredMort()], [setAlgaeParams()], [setDetritusParams()]
#'
#' @return  A MizerParams object with updated refuge parameters
#' @concept refugeParams
#' @export
setRefuge <- function(params, method, method_params = NULL,
                      # Parameters specific to each group
                      refuge_user = NULL, blocked_pred = NULL,
                      satiation = NULL,
                      # Parameters used by all methods
                      a_bar = NULL, b_bar = NULL,
                      w_settle = NULL, max_protect = NULL, tau = NULL,
                      use_dummy_fish_bins = TRUE,
                      ...) {
    # Object validation
    assert_that(is(params, "MizerParams"))

    # Find number of species for checks
    no_sp <- nrow(params@species_params)

    # === LENGTH-WEIGHT PARAMETER SETUP ===
    # Set default a_bar and b_bar values first
    if (is.null(a_bar)) {
        a_bar <- 0.025
    } else {
        if (!is.numeric(a_bar)) {
            stop("a_bar should be numeric.")
        }
        if (a_bar < 0) {
            stop("a_bar must be non-negative.")
        }
    }

    if (is.null(b_bar)) {
        b_bar <- 3
    } else {
        if (!is.numeric(b_bar)) {
            stop("b_bar should be numeric.")
        }
        if (b_bar < 0) {
            stop("b_bar must be non-negative.")
        }
    }

    # Check and create a and b columns if missing, set defaults for NAs
    missing_cols <- !c("a", "b") %in% names(params@species_params)
    if (any(missing_cols)) {
        if (missing_cols[1]) { # 'a' column missing
            params@species_params$a <- rep(NA_real_, nrow(params@species_params))
        }
        if (missing_cols[2]) { # 'b' column missing
            params@species_params$b <- rep(NA_real_, nrow(params@species_params))
        }
    }

    # Set defaults for missing a and b parameters using a_bar and b_bar
    if (anyNA(params@species_params[["a"]]) || anyNA(params@species_params[["b"]])) {
        missing_a <- is.na(params@species_params[["a"]])
        missing_b <- is.na(params@species_params[["b"]])

        if (any(missing_a)) {
            params@species_params[["a"]][missing_a] <- a_bar
        }
        if (any(missing_b)) {
            params@species_params[["b"]][missing_b] <- b_bar
        }

        if (any(missing_a) || any(missing_b)) {
            warning(
                "Missing values in species_params columns 'a' and/or 'b' have been set to average values (a_bar = ",
                a_bar, ", b_bar = ", b_bar, "). Consider providing species-specific length-weight parameters."
            )
        }
    }

    # === SPECIES PARAMETER CHECKS ===
    # Check that refuge_user is logical and the right length
    if (!("refuge_user" %in% colnames(params@species_params))) {
        if (is.null(refuge_user)) {
            warning("You have not provided values for refuge_user, so no species use refuge.")
            refuge_user <- rep(FALSE, no_sp)
        } else if (!is.logical(refuge_user)) {
            stop("The refuge_user values should be logical.")
        }
        if (length(refuge_user) != no_sp) {
            stop("refuge_user should have a value for every group.")
        }
        params@species_params$refuge_user <- refuge_user
    }

    # Check that blocked_pred is logical and the right length
    if (!("blocked_pred" %in% colnames(params@species_params))) {
        if (is.null(blocked_pred)) {
            warning("You have not provided values for blocked_pred, so all predators can access prey within refuge.")
            blocked_pred <- rep(FALSE, no_sp)
        } else if (!is.logical(blocked_pred)) {
            stop("The blocked_pred values should be logical.")
        }
        if (length(blocked_pred) != no_sp) {
            stop("blocked_pred should have a value for every group.")
        }
        params@species_params$blocked_pred <- blocked_pred
    }

    # Check that satiation is logical and the right length
    if (!("satiation" %in% colnames(params@species_params))) {
        if (is.null(satiation)) {
            # Calculate default satiation values based on feeding behavior
            # FALSE for carnivores (species that eat other species)
            # TRUE for resource consumers (only eat resources, not other species)

            # Check if species eat other species (row sums in interaction matrix > 0)
            interaction_rowsums <- rowSums(params@interaction)
            eats_other_species <- interaction_rowsums > 0

            # Check if species eat any resources
            resource_cols <- c(
                "resource_interaction", "interaction_algae", "interaction_detritus", "interaction_sponge",
                "interaction_encrusting", "interaction_massive", "interaction_cryptic", "interaction_branching"
            )

            existing_resource_cols <- resource_cols[resource_cols %in% names(params@species_params)]
            if (length(existing_resource_cols) > 0) {
                # Sum all resource interactions for each species (species are rows)
                resource_interaction_sums <- rowSums(params@species_params[, existing_resource_cols, drop = FALSE], na.rm = TRUE)
                eats_resources <- resource_interaction_sums > 0
            } else {
                eats_resources <- rep(FALSE, no_sp) # No resource interaction columns means no resource consumption
            }

            # Set satiation: FALSE for carnivores, TRUE for pure resource consumers
            satiation <- !eats_other_species & eats_resources

            warning("You have not specified whether species should have a satiation response. The default is FALSE for carnivores and TRUE for resource consumers.")
        } else if (!is.logical(satiation)) {
            stop("The satiation values should be logical.")
        }
        if (length(satiation) != no_sp) {
            stop("satiation should have a value for every group.")
        }
        params@species_params$satiation <- satiation
    }

    # === REFUGE METHOD VALIDATION ===
    # Check if the user provided one of the available refuge profile methods
    method_options <- c("sigmoidal", "binned", "competitive", "noncomplex")
    if (is.null(method)) {
        stop("You must provide the method to calculate the refuge profile.")
    } else if (!is.element(method, method_options)) {
        stop("Method must be 'sigmoidal','binned', 'competitive', 'noncomplex'.")
    }

    # === REFUGE PARAMETER SETUP ===
    # Set default values for parameters used by all methods

    # Minimum weight of fish protected by refuges at measured scale
    if (is.null(w_settle)) {
        w_settle <- 0.1
    } else {
        if (!is.numeric(w_settle)) {
            stop("w_settle should be numeric.")
        }
        if (w_settle < 0) {
            stop("w_settle must be non-negative.")
        }
    }

    # Maximum proportion of fish protected by refuge
    if (is.null(max_protect)) {
        max_protect <- 0.98
    } else {
        if (!is.numeric(max_protect)) {
            stop("max_protect should be numeric.")
        }
        if (max_protect < 0 || max_protect > 1) {
            stop("max_protect should be a proportion between 0 and 1")
        }
    }

    # Proportion of fish with access to refuge that are expected
    # to utilize it
    if (is.null(tau)) {
        tau <- 1
    } else {
        if (!is.numeric(tau)) {
            stop("tau should be numeric.")
        }
        if (tau < 0 || tau > 1) {
            stop("tau should be a proportion between 0 and 1")
        }
    }

    # Store all values directly in refuge_params slot
    params@refuge_params$method <- method
    params@refuge_params$a_bar <- a_bar
    params@refuge_params$b_bar <- b_bar
    params@refuge_params$w_settle <- w_settle
    params@refuge_params$max_protect <- max_protect
    params@refuge_params$tau <- tau
    params@refuge_params$use_dummy_fish_bins <- use_dummy_fish_bins

    #  method_params set up and checks
    if (method != "noncomplex") {
        # check if method_params provided
        if (is.null(method_params)) {
            stop("You must provide method specific parameters.")
        }

        # Convert list to data frame if needed
        if (is.list(method_params) && !is.data.frame(method_params)) {
            method_params <- as.data.frame(method_params)
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
        cnames <- colnames(method_params)

        ## Sigmoidal method
        if (method == "sigmoidal") {
            # Prop protect
            if (!("prop_protect" %in% cnames)) {
                stop("The sigmoidal method parameters dataframe needs a
                column called 'prop_protect' with the maximum proportion of
                fish to be protected.")
            } else if (method_params$prop_protect < 0 ||
                method_params$prop_protect > 1) {
                stop("prop_protect should be a proportion between 0 and 1.")
            }
            # L_refuge
            if (!("L_refuge" %in% cnames)) {
                stop("The sigmoidal method parameters dataframe needs a column
                called 'L_refuge' with the threshhold length (cm) for protected
                fish.")
            }
            # Slope
            if (is.null(method_params$slope)) {
                method_params$slope <- 100
            }
        }

        # Binned method
        if (method == "binned") {
            if (!("start_L" %in% cnames)) {
                stop("The binned method parameters dataframe needs a column called
                     'start_L' with the starting lengths (cm) for each size bin.")
            }
            if (!("end_L" %in% cnames)) {
                stop("The binned method parameters dataframe needs a column called
                     'end_L' with the end lengths (cm) for each size bin.")
            }
            if (!("prop_protect" %in% cnames)) {
                stop("The binned method parameters dataframe needs a column called
                     'prop_protect' with the proportion of fish protected
                     for each length bin.")
            } else if (any(method_params$prop_protect < 0) ||
                any(method_params$prop_protect > 1)) {
                stop("prop_protect should be a proportion between 0 and 1")
            }
            if (!all(method_params$start_L < method_params$end_L)) {
                stop("All bin start lengths must be less than bin end lengths.")
            }
        }

        # Competitive method
        if (method == "competitive") {
            if (!("start_L" %in% cnames)) {
                stop("The competitive method parameters dataframe needs a
                column called 'start_L' with the starting lengths (cm) for
                each size bin.")
            }
            if (!("end_L" %in% cnames)) {
                stop("The competitive method parameters dataframe needs a
                column called 'end_L' with the end lengths (cm) for
                each size bin.")
            }
            if (!("refuge_density" %in% cnames)) {
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


    # Store method_params in params object
    params@refuge_params[["method_params"]] <- as.data.frame(method_params)

    params@time_modified <- lubridate::now()

    return(params)
}

## Finds the refuge length bins by species and stores them in params
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
#' stored in a data frame called `refuge_lengths` in the `refuge_params` slot
#' of the `params` object.
#'
#' @inheritSection setRefuge Setting the refuge profile
#'
#' @param params a mizer params object
#' @param use_dummy_fish_bins   Logical. If FALSE, refuge thresholds and bin
#'                              boundaries are calculated using each species'
#'                              own length-weight parameters (`a`, `b`). If
#'                              TRUE (default), dummy fish parameters
#'                              (`a_bar`, `b_bar`) are used for conversions.
#'                              Applies to both sigmoidal and binned methods.
#'                              If not provided, uses value from
#'                              `params@refuge_params$species_specific_bins`.
#'                              See [setRefuge()] for details.
#'
#' @param ... Unused
#'
#' @return A mizer params object with updated refuge profiles
#' @concept refugeParams
#' @seealso [setRefuge()], [reefVulnerable()]
#' @export
getRefuge <- function(params, use_dummy_fish_bins = TRUE, ...) {
    # object - Check if mizerParams is valid
    assert_that(is(params, "MizerParams"))

    # Extract relevant data from params
    refuge_params <- params@refuge_params
    method_params <- params@refuge_params[["method_params"]]

    # Use value from params if not explicitly provided
    if (is.null(use_dummy_fish_bins)) {
        use_dummy_fish_bins <- isTRUE(refuge_params$use_dummy_fish_bins)
    }

    # Pull values from params
    w <- params@w
    sp <- params@species_params
    no_w <- length(params@w)
    no_sp <- dim(params@interaction)[1]

    # Set parameters used with all methods
    a_bar <- refuge_params$a_bar
    b_bar <- refuge_params$b_bar
    w_settle <- refuge_params$w_settle
    max_protect <- refuge_params$max_protect

    # Store which functional groups use refuge
    refuge_user <- sp$refuge_user

    # Noncomplex reef ----------------------------------------------------------
    if (refuge_params$method == "noncomplex") {
        # Create matrix to store proportions for each species
        refuge <- matrix(0, nrow = no_sp, ncol = no_w)
        rownames(refuge) <- rownames(params@initial_n)
        colnames(refuge) <- colnames(params@initial_n)
        # store refuge and bin indices in params object
        params@refuge_params$refuge <- refuge

        # Sigmoidal method ---------------------------------------------------------
    } else if (refuge_params$method == "sigmoidal") {
        # Pull slope and proportion of fish to be protected from method_params
        prop_protect <- method_params$prop_protect
        slope <- method_params$slope
        # Create matrix to store proportions for each species
        refuge <- matrix(0, nrow = no_sp, ncol = no_w)
        rownames(refuge) <- rownames(params@initial_n)
        colnames(refuge) <- colnames(params@initial_n)
        if (isFALSE(use_dummy_fish_bins)) {
            # Species-specific threshold: use each species' a and b
            for (i in 1:no_sp) {
                W_refuge_i <- sp[["a"]][i] * method_params$L_refuge^sp[["b"]][i]
                denom <- 1 + exp(slope * (w - W_refuge_i))
                ref <- ifelse(w > w_settle, prop_protect / denom, 0)
                ref[ref > max_protect] <- max_protect
                refuge[i, ] <- refuge_user[i] * ref
            }
            refuge_lengths <- data.frame(species = sp$species, L_refuge = method_params$L_refuge)
        } else if (isTRUE(use_dummy_fish_bins)) {
            # Dummy fish threshold: use a_bar and b_bar
            W_refuge <- a_bar * method_params$L_refuge^b_bar
            denom <- 1 + exp(slope * (w - W_refuge))
            ref <- ifelse(w > w_settle, prop_protect / denom, 0)
            ref[ref > max_protect] <- max_protect
            for (i in 1:no_sp) {
                refuge[i, ] <- refuge_user[i] * ref
            }
            L_refuge.i <- (W_refuge / sp[["a"]])^(1 / sp[["b"]])
            refuge_lengths <- data.frame(species = sp$species, L_refuge = L_refuge.i)
        }
        params@refuge_params$refuge <- refuge
        params@refuge_params$refuge_lengths <- refuge_lengths
        params@time_modified <- lubridate::now()

        # Binned method ------------------------------------------------------------
    } else if (refuge_params$method == "binned") {
        # Initialize storage
        ref <- rep(0, no_w)
        start_l.i <- list(1)
        end_l.i <- list(1)
        bin.id <- list(1)
        no_bins <- nrow(method_params)
        if (isFALSE(use_dummy_fish_bins)) {
            # Species-specific bin boundaries
            for (k in 1:no_bins) {
                for (i in 1:no_sp) {
                    start_w <- sp[["a"]][i] * method_params$start_L[[k]]^sp[["b"]][i]
                    end_w <- sp[["a"]][i] * method_params$end_L[[k]]^sp[["b"]][i]
                    start_w[start_w < w_settle] <- w_settle
                    bin.id[[paste0("sp", i, "_bin", k)]] <- which(w >= start_w & w <= end_w)
                    # Assign prop_protect for this bin/species
                    ref[bin.id[[paste0("sp", i, "_bin", k)]]] <- method_params$prop_protect[k]
                    # Calculate length bins for each species
                    start_l.i[[paste0("sp", i, "_bin", k)]] <- (start_w / sp[["a"]][i])^(1 / sp[["b"]][i])
                    end_l.i[[paste0("sp", i, "_bin", k)]] <- (end_w / sp[["a"]][i])^(1 / sp[["b"]][i])
                }
            }
        } else {
            # Dummy fish bin boundaries
            for (k in 1:no_bins) {
                start_w <- a_bar * method_params$start_L[[k]]^b_bar
                end_w <- a_bar * method_params$end_L[[k]]^b_bar
                start_w[start_w < w_settle] <- w_settle
                bin.id[[k]] <- which(w >= start_w & w <= end_w)
                ref[bin.id[[k]]] <- method_params$prop_protect[k]
                # Calculate length bins for each species
                start_l.i[[k]] <- (start_w / sp[["a"]])^(1 / sp[["b"]])
                names(start_l.i)[[k]] <- c(paste("start", k, sep = ""))
                end_l.i[[k]] <- (end_w / sp[["a"]])^(1 / sp[["b"]])
                names(end_l.i)[[k]] <- c(paste("end", k, sep = ""))
            }
        }
        # Create matrix to store proportions for each species
        refuge <- matrix(rep(ref), nrow = no_sp, ncol = no_w, byrow = TRUE)
        rownames(refuge) <- rownames(params@initial_n)
        colnames(refuge) <- colnames(params@initial_n)
        # Make sure none of the values are higher than maximum protection allowed
        refuge[refuge > max_protect] <- max_protect
        # Account for species that don't utilize refuge
        refuge <- refuge_user * refuge
        # store refuge and bin indices in params object
        params@refuge_params$refuge <- refuge
        params@refuge_params$bin.id <- bin.id
        # store length bins by functional group in params object
        start_l.i <- t(do.call(rbind, start_l.i))
        end_l.i <- t(do.call(rbind, end_l.i))
        refuge_lengths <- cbind(start_l.i, end_l.i)
        row.names(refuge_lengths) <- sp$species
        params@refuge_params$refuge_lengths <- refuge_lengths
        params@time_modified <- lubridate::now()

        ## Competitive method ------------------------------------------------------
    } else if (refuge_params$method == "competitive") {
        # Initialize empty list to hold number of competitors for each bin
        bin.id <- list(1)
        start_l.i <- list(1)
        end_l.i <- list(1)
        no_bins <- nrow(method_params)
        if (isFALSE(use_dummy_fish_bins)) {
            # Use species-specific a/b to convert length bins to weights
            for (k in 1:no_bins) {
                for (i in 1:no_sp) {
                    start_w <- sp[["a"]][i] * method_params$start_L[[k]]^sp[["b"]][i]
                    end_w <- sp[["a"]][i] * method_params$end_L[[k]]^sp[["b"]][i]
                    start_w[start_w < w_settle] <- w_settle
                    bin.id[[paste0("sp", i, "_bin", k)]] <- which(params@w >= start_w & params@w <= end_w)
                    # Calculate length bins for each species
                    start_l.i[[paste0("sp", i, "_bin", k)]] <- (start_w / sp[["a"]][i])^(1 / sp[["b"]][i])
                    end_l.i[[paste0("sp", i, "_bin", k)]] <- (end_w / sp[["a"]][i])^(1 / sp[["b"]][i])
                }
            }
        } else {
            # Use dummy fish parameters to convert length bins to weights, then convert to species-specific lengths
            for (k in 1:no_bins) {
                start_w <- a_bar * method_params$start_L[[k]]^b_bar
                end_w <- a_bar * method_params$end_L[[k]]^b_bar
                start_w[start_w < w_settle] <- w_settle
                bin.id[[k]] <- which(params@w >= start_w & params@w <= end_w)
                # Calculate length bins for each species
                start_l.i[[k]] <- (start_w / sp[["a"]])^(1 / sp[["b"]])
                names(start_l.i)[[k]] <- c(paste("start", k, sep = ""))
                end_l.i[[k]] <- (end_w / sp[["a"]])^(1 / sp[["b"]])
                names(end_l.i)[[k]] <- c(paste("end", k, sep = ""))
            }
        }
        # Store indices of each bin
        params@refuge_params$bin.id <- bin.id
        # Store length bins by species group in a data frame
        start_l.i <- t(do.call(rbind, start_l.i))
        end_l.i <- t(do.call(rbind, end_l.i))
        refuge_lengths <- cbind(start_l.i, end_l.i)
        row.names(refuge_lengths) <- sp$species
        params@refuge_params$refuge_lengths <- refuge_lengths
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
#' @param new_refuge    A boolean value that states whether this new refuge
#'                      profile is being used for simulation. Determines
#'                      whether algae and detritus production are tuned when
#'                      the model is run to steady state. Defaults to FALSE.
#'
#' @param new_method    The new method to be used for setting the refuge
#'                      profile. Options are "sigmoidal", "binned",
#'                      "competitive", or "noncomplex". If no method is
#'                      provided, this defaults to the same method as is
#'                      currently being used in the simulation.
#'
#' @param new_method_params  A data frame or list specifying parameters required for
#'                           the chosen method. For 'sigmoidal', must include 'L_refuge'
#'                           and 'prop_protect'. For 'binned' and 'competitive', must
#'                           include 'start_L', 'end_L', and 'prop_protect'
#'                           (or 'refuge_density' for competitive). Only necessary if
#'                           changing methods.
#' @param scale_bin To be used with the "binned" or "competitive" method only. A number
#'                  or vector (length = number of bins) to multiply the `prop_protect`
#'                  (binned) or `refuge_density` (competitive) values for each size bin.
#'                  Changes the proportion of fish protected. If a single value is given,
#'                  it is applied to all bins.
#' @details At least one of the arguments new_method, new_method_params, new_L_refuge,
#'          new_prop_protect, or scale_bin must be provided. For 'binned' and 'competitive'
#'          methods, scale_bin must be a value or vector matching the number of bins.
#'
#' @param new_L_refuge  To be used with "sigmoidal" method only. The new value
#'                      for the length at which refuge becomes scarce in cm.
#'
#' @param new_prop_protect  To be used with "sigmoidal" method only. The new
#'                          value for the maximum proportion of fish to protect.
#'
#' @param ... Unused
#'
#' @return A mizer object with updated refuge profiles
#' @concept refugeParams
#' @export
newRefuge <- function(params, new_refuge = FALSE,
                      # Fully changing method
                      new_method = NULL, new_method_params = NULL,
                      # Sigmoidal - changing refuge length prop protect
                      new_L_refuge = NULL, new_prop_protect = NULL,
                      # Binned - scaling bin prop
                      scale_bin = NULL, ...) {
    # Check that the user provided at least one new input
    inputs <- list(
        new_method, new_method_params,
        new_L_refuge, new_prop_protect, scale_bin
    )
    if (all(sapply(inputs, is.null))) {
        stop("Error: At least one input must be provided.")
    }

    # Save new refuge flag
    params@other_params$new_refuge <- new_refuge

    # object check
    # Check if given mizerParams object is valid
    assert_that(is(params, "MizerParams"))

    # Extract relevant data from params
    refuge_params <- params@refuge_params
    method_params <- params@refuge_params[["method_params"]]

    # method checks
    # Check if the user provided one of the available refuge profile methods
    m_opt <- c("sigmoidal", "binned", "competitive", "noncomplex")
    if (is.null(new_method)) {
        warning("Since you did not provide a new method to set the refuge profile,
                    I will use the one currently stored in params@refuge_params$method.")
        # If user did not provide a method, use old one
        new_method <- refuge_params$method
        # If user did provide a method, check that it's one of the options
    } else if (!is.element(new_method, m_opt)) {
        stop("Method must be 'sigmoidal','binned', 'competitive', 'noncomplex'.")
    }

    #  method_params checks
    if (new_method != "noncomplex") {
        # check if method_params provided
        if (is.null(new_method_params)) {
            # Return an error if they are trying to switch to methods
            if (refuge_params$method != new_method) {
                stop("You must provide a new method_params data frame to
                         switch methods.")
            }
            # Check method
            ## Sigmoidal
            if (new_method == "sigmoidal") {
                # Check if new L or new_prop given
                if (all(is.null(new_L_refuge), is.null(new_prop_protect))) {
                    stop("You must provide either a new L_refuge, a new
                             prop_protect, or a new method_params data frame.")
                }
                # If new L_refuge given, check that it's positive
                if (!is.null(new_L_refuge)) {
                    if (new_L_refuge < 0) {
                        stop("new_L_refuge must be non-negative.")
                    }
                    # If it's not given, use the old one
                } else {
                    new_L_refuge <- method_params$L_refuge
                }
                # If new prop_protect given, check that it's between 0 and 1
                if (!is.null(new_prop_protect)) {
                    if (method_params$prop_protect < 0 ||
                        method_params$prop_protect > 1) {
                        stop("prop_protect should be a proportion between 0 and 1.")
                    }
                    # If it's not given, use the old one
                } else {
                    
                    new_prop_protect <- method_params$prop_protect
                }

                new_mp <- method_params
                new_mp$L_refuge <- new_L_refuge
                new_mp$prop_protect <- new_prop_protect

                ## Binned or Competitive
            } else if (new_method == "binned" || new_method == "competitive") {
                # Find number of bins used in old method
                no_bins <- length(method_params)
                if (is.null(scale_bin)) {
                    stop("You must provide either a value or vector of values
                             for scale_bin or a new method_params data frame.")
                }
                # Make sure scaling vector is the right length
                if (!(length(scale_bin) == 1) ||
                    !(length(scale_bin) == no_bins)) {
                    stop("scale_bin must have length 1 or have an value
                             for every bin.")
                }
                # Check that scale_bin_prop is positive
                if (scale_bin < 0) {
                    stop("scale_bin must be non-negative.")
                }
                # Calculate new bins
                if (new_method == "binned") {
                    new_mp <- method_params
                    new_mp$prop_protect <- scale_bin * new_mp$prop_protect
                } else {
                    new_mp <- method_params
                    new_mp$refuge_density <- scale_bin * new_mp$refuge_density
                }
            }
        } else {
            if (is.list(new_method_params) && !is.data.frame(new_method_params)) {
                new_method_params <- as.data.frame(new_method_params)
            }
            new_mp <- new_method_params
        }
    }

    # Update parameters
    # Store new parameters
    params <- setRefuge(
        params = params,
        method = new_method,
        method_params = new_mp
    )

    # Find new refuge profile
    params <- getRefuge(params = params)

    # Save time parameters were modified
    params@time_modified <- lubridate::now()

    # Return mizer params object
    return(params)
}


#' Prepare a steady state model for projections with degradation
#'
#' This function stores degradation parameters in the `refuge_params` slot of
#' the params object, and algae adjustment parameters in the `algae_params` slot.
#'
#' @param params a mizer object
#'
#' @param trajectory    Optional character string naming the degradation scenario.
#'                      Used for documentation and identification. No functional
#'                      effect on the simulation.
#'
#' @param deg_scale     A matrix (refuge bins x time) giving scaling factors for
#'                      refuge density. Column 1 represents the bleaching year
#'                      (initial impact), and subsequent columns represent years
#'                      1, 2, 3... post-bleaching. Values are multiplied by the
#'                      previous timestep's refuge density. Default scaling matrices
#'                      for 15 years with "rubble", "algae", and "recovery"
#'                      trajectories are included as data objects in the package.
#'
#' @param bleach_time   The year of the simulation to implement bleaching.
#'                      Defaults to year 2.
#'
#' @param algae_boost  Logical. Should algae growth and/or carrying capacity
#'                      be adjusted in response to bleaching? Default is FALSE.
#'
#' @param algae_growth_boost Numeric vector. Multipliers for algae growth rate
#'                      at bleaching (element 1) and during post-bleaching years
#'                      (subsequent elements). Only used if algae_boost = TRUE.
#'                      Vector length determines duration: length 4 means bleaching
#'                      year plus 3 post-bleaching years. Default is
#'                      c(1.11, 1.11, 1.11, 1.11) for 4 years of boosting.
#'
#' @param algae_capacity_boost Numeric vector. Multipliers for algae carrying
#'                      capacity at bleaching (element 1) and during post-bleaching
#'                      years (subsequent elements). Only used if algae_boost = TRUE.
#'                      If shorter than algae_growth_boost, will be padded with 1s
#'                      (no change). Default is c(2.0) (only boosts at bleaching year).
#'
#' @param ... Unused
#'
#' @return A mizer object with updated degradation parameters stored in `refuge_params`
#'         and algae adjustment parameters in `algae_params`.
#' @concept degradation
#' @seealso [algae_dynamics_cc()],[detritus_dynamics_cc()],
#'          [tuneUR_cc()], [reefDegrade()]
#' @export
setDegradation <- function(params, trajectory = NULL, deg_scale,
                           bleach_time = 2,
                           degrade = FALSE,
                           algae_boost = FALSE,
                           algae_growth_boost = c(1.11, 1.11, 1.11, 1.11),
                           algae_capacity_boost = c(2.0),
                           ...) {
    # Validate deg_scale
    if (!is.matrix(deg_scale) && !is.data.frame(deg_scale)) {
        stop("deg_scale must be a matrix or data frame.")
    }
    deg_scale <- as.matrix(deg_scale)

    if (ncol(deg_scale) < 1) {
        stop("deg_scale must have at least one column (bleaching year impact).")
    }

    # Validate algae parameters if adjustment is enabled
    if (algae_boost) {
        if (!is.numeric(algae_growth_boost) || any(algae_growth_boost <= 0)) {
            stop("algae_growth_boost must be a positive numeric vector.")
        }
        if (!is.numeric(algae_capacity_boost) || any(algae_capacity_boost <= 0)) {
            stop("algae_capacity_boost must be a positive numeric vector.")
        }
        # Pad capacity boost with 1s (no change) if shorter than growth boost
        if (length(algae_capacity_boost) < length(algae_growth_boost)) {
            algae_capacity_boost <- c(
                algae_capacity_boost,
                rep(1, length(algae_growth_boost) - length(algae_capacity_boost))
            )
        }
    }

    # Add new parameters to params object
    params@refuge_params$degrade <- degrade
    params@refuge_params$t_bleach <- bleach_time
    params@refuge_params$trajectory <- trajectory
    params@refuge_params$deg_scale <- deg_scale

    # Store algae adjustment parameters
    params@algae_params$algae_boost <- algae_boost
    if (algae_boost) {
        params@algae_params$algae_growth_boost <- algae_growth_boost
        params@algae_params$algae_capacity_boost <- algae_capacity_boost
    }

    # Save time parameters were modified
    params@time_modified <- lubridate::now()

    return(params)
}

#' Switch to unstructured resource dynamics with carrying capacities
#'
#' This function calculates new carrying capacities for algae and detritus
#' resources and switches the resource dynamics from the standard linear
#' dynamics to fluxes with a carrying capacity. The carrying capacity is set
#' at two times the current steady state biomass of the resource.
#'
#' See [algae_dynamics_cc()] and [detritus_dynamics_cc()] for additional detail.
#'
#' @param params a \linkS4class{MizerParams} object
#' @param cap   Numeric. Value to scale the steady state biomass by. Defaults to
#'              1.5, setting the carrying capacity 50% higher than the current
#'              standing biomass.
#' @param ... Unused
#'
#' @return A mizer object with updated degradation parameters profiles
#' @concept Uresources
#' @seealso [algae_dynamics_cc()],[detritus_dynamics_cc()],
#'          [tuneUR_cc()]
#' @export
setURcapacity <- function(params, cap = 1.5, ...) {
    # Change capacity boolean to true
    params@other_params$use_UR_cc <- TRUE

    # Calculate algae and detritus carrying capacities
    # from steady state biomasses
    # Algae
    ba <- algae_biomass(params)
    new_a_carry <- cap * ba
    params@algae_params$algae_capacity <- new_a_carry
    # Detritus
    bd <- detritus_biomass(params)
    new_d_carry <- cap * bd
    params@detritus_params$detritus_capacity <- new_d_carry

    # Switch to carrying capacity dynamics
    params@other_dynamics$algae <- "algae_dynamics_cc"
    params@other_dynamics$detritus <- "detritus_dynamics_cc"

    # Save time parameters were modified
    params@time_modified <- lubridate::now()

    return(params)
}

# These functions will be useful when mizerReef is ready to project into the
# future. We are not there yet.

# A Version of getEGrowth that allows for mizerSim/ gets values through time ---
# Get energy rate available for growth at ALL time steps
#
# Calculates the energy rate \eqn{g_i(w)} (grams/year) available for growth
# by functional group and size after metabolism, movement and reproduction
# have been accounted for.
#
# This function replace's mizer's original getEgrowth
#
# @inherit mizerEGrowth
# @inheritParams mizerRates
#
# @return A two dimensional array (prey species x prey size)
# @export
# @seealso [getERepro()], [getEReproAndGrowth()]
# @family rate functions
# @examples
# \donttest{
# params <- NS_params
# # Project with constant fishing effort for all gears for 20 time steps
# sim <- project(params, t_max = 20, effort = 0.5)
# # Get the energy at a particular time step
# getEGrowth(params, n = N(sim)[15, , ], n_pp = NResource(sim)[15, ], t = 15)
# }
# getEGrowth <- function(params, n = initialN(params),
#                        n_pp = initialNResource(params),
#                        n_other = initialNOther(params),
#                        t = 0, ...) {
#     params <- validParams(params)
#     f <- get(params@rates_funcs$EGrowth)
#     g <- f(params, n = n, n_pp = n_pp, n_other = n_other, t = t,
#            e_repro = getERepro(params, n = n, n_pp = n_pp,
#                                n_other = n_other, t = t),
#            e = getEReproAndGrowth(params, n = n, n_pp = n_pp,
#                                   n_other = n_other, t = t))
#     dimnames(g) <- dimnames(params@metab)
#     g

# Productivity through time ----------------------------------------------------
# getProductivity <- function(object,
#                             min_fishing_l = 7,...) {
#     if (is(object, "MizerSim")) {
#         sim <- object
#         size_range <- get_size_range_array(sim@params,
#                                            min_l = min_fishing_l...)
#         g <- getEGrowth(object)
#         productivity <- params@initial_n * size_range
#         return(
#     }

# # Get refuge profile through time --------------------------------------------
# @param object A `MizerParams` object or a `MizerSim` object
# @inheritParams reefRates
# @inheritParams get_time_elements
#
# @return If a `MizerParams` object is passed in, the function returns a two
#   dimensional array (prey species x prey size) based on the refuge parameters
#   and abundances also passed in.
#
#   If  a `MizerSim` object is passed in, the function returns a three
#   dimensional array (time step x prey species x prey size) with the
#   vulnerability calculated at every time step in the simulation.
#
#   Currently, both abundance and diversity of refuges stay static through
#   time. The simple and binned methods stay constant while the data method
#   is density dependent.
# getRefuge <- function(object,
#                           n, n_pp, n_other,
#                           time_range, drop = FALSE, ...) {
#     if (is(object, "MizerParams")) {
#
#         params <- validParams(object)
#         if (missing(time_range)) time_range <- 0
#         t <- min(time_range)
#         if (missing(n)) n <- params@initial_n
#         if (missing(n_pp)) n_pp <- params@initial_n_pp
#         if (missing(n_other)) n_other <- params@initial_n_other
#
#         # Calculate vulnerability
#         f <- get(params@rates_funcs$Vulnerable)
#         vulnerable <- f(params, n = n, n_pp = n_pp, n_other = n_other, t = t)
#         dimnames(vulnerable) <- dimnames(params@metab)
#         return(vulnerable)
#
#     } else {
#
#         sim <- object
#         if (missing(time_range)) {
#             time_range <- dimnames(sim@n)$time
#         }
#
#         time_elements <- get_time_elements(sim, time_range)
#
#         vul_time <- plyr::aaply(which(time_elements), 1, function(x) {
#             # Necessary as we only want single time step but may only have 1
#             # species which makes using drop impossible
#             n <- array(sim@n[x, , ], dim = dim(sim@n)[2:3])
#             dimnames(n) <- dimnames(sim@n)[2:3]
#             n_other <- sim@n_other[x, ]
#             names(n_other) <- dimnames(sim@n_other)$component
#             t <- as.numeric(dimnames(sim@n)$time[[x]])
#             vul <- getVulnerable(sim@params, n = n,
#                                  n_pp = sim@n_pp[x, ],
#                                  n_other = n_other,
#                                  time_range = t)
#             return(vul)
#         }, .drop = FALSE)
#
#         # Before we drop dimensions we want to set the time dimname
#         names(dimnames(vul_time))[[1]] <- "time"
#         vul_time <- vul_time[, , , drop = drop]
#         return(vul_time)
#
#     }
# }

# Plot vulnerability through time ----------------------------------------------

# Plot producitivity through time ----------------------------------------------

# Calculate initial vulnerability to predation
#
# This function uses the model parameters and refuge parameters to calculate
# initial values for the vulnerability matrix and stores it in the
# `other_params` slot of the `params` object.
#
# Must be run after the `setRefuge()` function.
#
# @param params The model parameters. An object of type
#   \linkS4class{MizerParams}.
# @param ... not used
#
# @export
# @concept helper
# @return A two dimensional array (prey species x prey size) listing the
# proportion of fish at that size that are vulnerable to being encountered by
# predators.
# get_initial_vulnerable <- function(params, ...) {
#
#     # Extract relevant data from params
#     refuge_params <- params@other_params[['refuge_params']]
#     method_params <- params@other_params[['method_params']]
#
#     # Pull values from params
#     w <- params@w
#     no_w <- length(params@w)
#     no_sp <- dim(params@interaction)[1]
#
#     # Initialize storage for the array of refuge proportions
#     refuge <- matrix(0, nrow = no_sp, ncol = no_w)
#     rownames(refuge) <- rownames(params@initial_n)
#     colnames(refuge) <- colnames(params@initial_n)
#
#     # Set parameters used with all methods
#     min_ref_w   <- params@species_params$min_ref_w
#     max_protect <- refuge_params$max_protect
#     tau         <- refuge_params$tau
#
#     # Store which functional groups use refuge
#     refuge_user <- params@species_params$refuge_user
#
#     ### SIMPLE METHOD ---------------------------------------------------------
#     if (refuge_params$method == "simple"){
#
#         # Pull slop and proportion of fish to be protected from method_params
#         prop_protect <- method_params$prop_protect
#         slope <- method_params$slope
#
#         # Convert length to weight
#         max_W <- params@species_params[["a"]] *
#             method_params$max_L ^ params@species_params[["b"]]
#
#         # Find indices of fish in size range to protect
#         # Set threshold weight - no organisms smaller than min_ref_length
#         # can utilize refuge to escape predators
#         # TRUE for fish larger than minimum protected size
#         min <- t(sapply(min_ref_w, function(x) params@w >= x))
#         # TRUE for fish smaller than maximum protected weight
#         max <- t(sapply(max_W, function(x) params@w <= x))
#         # TRUE for fish that meet both conditions
#         bin_fish <- min + max
#         idx.bin <- which(bin_fish == TRUE)
#
#         # Calculate protection level for fish
#         refuge[idx.bin] <- (-1 * prop_protect) /
#             (1 + exp(-1 * slope*(w - max_W))) + prop_protect
#
#         ### BINNED METHOD ----------------------------------------------------------
#     } else if (refuge_params$method == "binned") {
#
#         # Loop through each refuge bin
#         for (k in 1:nrow(method_params)) {
#
#             # Calculate start and end of weight bins for each functional group
#             # based on unique as and bs
#             start_w <- params@species_params[["a"]] *
#                 method_params$start_L[[k]] ^ params@species_params[["b"]]
#
#             end_w  <- params@species_params[["a"]] *
#                 method_params$end_L[[k]] ^ params@species_params[["b"]]
#
#             # Set threshold weight - no organisms smaller than min_ref_length
#             # can utilize refuge to escape predators
#             idx.sm <- which(start_w < min_ref_w)
#             start_w[idx.sm] <- min_ref_w[idx.sm]
#
#             # Find indices of fish in size range to protect
#             # TRUE for fish larger than start weight of bin
#             start_n <- t(sapply(start_w, function(x) params@w >= x))
#             # TRUE for fish smaller than end weight of bin
#             end_n <- t(sapply(end_w, function(x) params@w <= x))
#             # TRUE for fish that meet both conditions
#             bin_fish <- start_n & end_n
#
#             # Find indices of functional group x weight that are in size bin k
#             idx.bin = which(bin_fish == TRUE)
#
#             # Set vulnerability for fish in size bin to provided value
#             refuge[idx.bin] = method_params$prop_protect[k]
#         }
#
#         ### DATA METHOD ------------------------------------------------------------
#     } else if (refuge_params$method == "data") {
#
#         # Initialize empty list to hold number of competitors for each bin
#         competitor_density = numeric(nrow(method_params))
#
#         # Loop through each refuge bin
#         for (k in 1:nrow(method_params)) {
#
#             # Calculate start and end of weight bins for each
#             # functional group based on unique as and bs
#             start_w <- params@species_params[["a"]] *
#                 method_params$start_L[[k]] ^ params@species_params[["b"]]
#
#             end_w  <- params@species_params[["a"]] *
#                 method_params$end_L[[k]] ^ params@species_params[["b"]]
#
#             # Set threshold weight - no organisms smaller than
#             # minimum size that can utilize refuge
#             idx.sm <- which(start_w < min_ref_w)
#             start_w[idx.sm] <- min_ref_w[idx.sm]
#
#             # Calculate competitor density - number of fish that use refuge in
#             # each size bin
#
#             # Create matrix of boolean values indicating whether fish is
#             # in size range of bin k
#             # TRUE for fish larger than start weight
#             start_n <- t(sapply(start_w, function(x) params@w >= x))
#             # TRUE for fish smaller than end weight
#             end_n <- t(sapply(end_w, function(x) params@w <= x))
#             # TRUE for fish that meet both conditions
#             bin_fish <- start_n & end_n
#
#             # Number of competitors from each functional group in bin
#             competitors <- (initial_n * bin_fish) %*% params@dw
#
#             # Eliminate functional groups that don't use refuge and sum
#             competitor_density[k] <- sum(refuge_user * competitors)
#
#             # Find indices of fish within size bin k
#             idx.bin = which(bin_fish == TRUE)
#
#             # Set vulnerability for fish in size bin based on the number of
#             # available refuges and the number of competitors
#             refuge[idx.bin] <- ifelse(competitor_density[k] == 0, max_protect,
#                                       tau * method_params$refuge_density[k] /
#                                           competitor_density[k])
#         }
#     }
#
#     # Make sure none of the values are higher than maximum protection allowed
#     refuge[refuge > max_protect] = max_protect
#
#     # Account for species that don't utilize refuge
#     vul = 1 - (refuge_user*refuge)
#
#     # Store in params
#     params@other_params[['initial_vulnerable']] <- as.data.frame(vul)
#
#     params@time_modified <- lubridate::now()
#
#     return(params)
# }
#
#
#

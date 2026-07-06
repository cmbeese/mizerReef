# New params
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
    params@other_params$refuge_params$method <- method
    params@other_params$refuge_params$a_bar <- a_bar
    params@other_params$refuge_params$b_bar <- b_bar
    params@other_params$refuge_params$w_settle <- w_settle
    params@other_params$refuge_params$max_protect <- max_protect
    params@other_params$refuge_params$tau <- tau
    params@other_params$refuge_params$use_dummy_fish_bins <- use_dummy_fish_bins

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
    params@other_params$refuge_params[["method_params"]] <- as.data.frame(method_params)

    params@time_modified <- lubridate::now()

    return(params)
}

# Old params
setRefuge <- function(params, method, method_params = NULL,
                      # Parameters specific to each group
                      refuge_user = NULL, bad_pred = NULL,
                      satiation = NULL,
                      # Parameters used by all methods
                      a_bar = NULL, b_bar = NULL,
                      w_settle = NULL, max_protect = NULL, tau = NULL, ...) {
    # object check
    # Check if given mizerParams object is valid
    assert_that(is(params, "MizerParams"))

    # Check that a and b parameters are present for all species -
    # needed for l2w conversion
    if (any(!c("a", "b") %in% names(params@species_params))) {
        stop(
            "species_params slot must have columns 'a' and 'b' for ",
            "length-weight conversion"
        )
    }
    if (anyNA(params@species_params[["a"]]) ||
        anyNA(params@species_params[["b"]])) {
        message("There are NAs in the species_params
                    columns 'a' and 'b'. You must provide these parameters
                    for these species if you want them to use refuge.")
    }

    # Find number of species for checks
    no_sp <- nrow(params@species_params)

    # species_params checks
    # Check that refuge_user is logical and the right length
    if (!("refuge_user" %in% colnames(params@species_params))) {
        if (is.null(refuge_user)) {
            stop("You need to provide values for refuge_user")
        } else if (!is.logical(refuge_user)) {
            stop("The refuge_user values should be logical.")
        }
        if (length(refuge_user) != no_sp) {
            stop("refuge_user should have a value for every group.")
        }
        params@species_params$refuge_user <- refuge_user
    }

    # Check that bad_pred is logical and the right length
    if (!("bad_pred" %in% colnames(params@species_params))) {
        if (is.null(bad_pred)) {
            stop("You need to provide values for bad_pred")
        } else if (!is.logical(bad_pred)) {
            stop("The bad_pred values should be logical.")
        }
        if (length(refuge_user) != no_sp) {
            stop("bad_pred should have a value for every group.")
        }
        params@species_params$bad_pred <- bad_pred
    }

    # Check that satiation is logical and the right length
    if (!("satiation" %in% colnames(params@species_params))) {
        if (is.null(satiation)) {
            stop("You need to provide values for satiation")
        } else if (!is.logical(satiation)) {
            stop("The satiation values should be logical.")
        }
        if (length(satiation) != no_sp) {
            stop("satiation should have a value for every group.")
        }
        params@species_params$satiation <- satiation
    }

    # refuge_params set up and checks
    # Check if the user provided one of the available refuge profile methods
    method_options <- c("sigmoidal", "binned", "competitive", "noncomplex")
    if (is.null(method)) {
        stop("You must provide the method to calculate the refuge profile.")
    } else if (!is.element(method, method_options)) {
        stop("Method must be 'sigmoidal','binned', 'competitive', 'noncomplex'.")
    }

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

    # a_bar for fish dummies
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

    # a_bar for fish dummies
    if (is.null(b_bar)) {
        b_bar <- 3
    } else {
        if (!is.numeric(a_bar)) {
            stop("b_bar should be numeric.")
        }
        if (a_bar < 0) {
            stop("b_bar must be non-negative.")
        }
    }

    # Store all values in refuge_params data frame
    refuge_params <- data.frame(
        method, a_bar, b_bar,
        w_settle, max_protect, tau
    )

    #  method_params set up and checks
    if (method != "noncomplex") {
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
        cnames <- colnames(method_params)

        ## Sigmoidal method
        if (refuge_params$method == "sigmoidal") {
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
        if (refuge_params$method == "binned") {
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
        if (refuge_params$method == "competitive") {
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


    # Store in params object --
    params@other_params[["refuge_params"]] <- as.data.frame(refuge_params)

    params@other_params[["method_params"]] <- as.data.frame(method_params)

    params@time_modified <- lubridate::now()

    return(params)
}

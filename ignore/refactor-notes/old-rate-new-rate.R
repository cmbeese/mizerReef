# New reef vulnerable for competitive method
reefVulnerable <- function(params, n, n_pp, n_other, t, new_rd = NULL, ...) {
    # Extract relevant data from params
    method_params <- params@refuge_params$method_params
    degrade <- params@refuge_params$degrade

    # If degradation is being implemented, calculate new refuge density
    if (degrade == TRUE) {
        if (is.null(new_rd)) {
            new_rd <- reefDegrade(params, n, n_pp, n_other, t, ...)
        }
    } else {
        new_rd <- method_params$refuge_density
    }

    # Set parameters used with all methods
    max_protect <- params@refuge_params$max_protect
    tau <- params@refuge_params$tau

    # Pull no of species and size bins
    no_w <- length(params@w)
    no_sp <- dim(params@interaction)[1]

    # Store which functional groups use refuge
    refuge_user <- params@species_params$refuge_user

    # Static methods -----------------------------------------------------------
    static <- c("sigmoidal", "binned", "noncomplex")

    if (params@refuge_params$method %in% static) {
        refuge <- params@refuge_params$refuge
        vulnerable <- 1 - refuge

        # Competitive method -------------------------------------------------------
    } else if (params@refuge_params$method == "competitive") {
        # Determine bin.id structure
        bin_id_list <- params@refuge_params$bin.id
        # Check if bin.id uses character keys (species-specific bins)
        bin_names <- names(bin_id_list)
        use_dummy_fish_bins <- !is.null(bin_names) && !any(grepl("sp", bin_names))

        # Initialize storage for the array of refuge proportions
        refuge <- matrix(0, nrow = no_sp, ncol = no_w)
        rownames(refuge) <- rownames(params@initial_n)
        colnames(refuge) <- colnames(params@initial_n)

        if (!use_dummy_fish_bins) {
            # Handle species-specific bins
            competitor_density <- numeric(length(bin_id_list))
            bin_keys <- bin_names
            for (idx in seq_along(bin_keys)) {
                key <- bin_keys[idx]
                # Parse species and bin from key
                sp_match <- regmatches(key, regexec("sp([0-9]+)_bin([0-9]+)", key))[[1]]
                if (length(sp_match) == 3) {
                    i <- as.integer(sp_match[2])
                    k <- as.integer(sp_match[3])
                    bin.id <- bin_id_list[[key]]
                    # Create logical vector and use to get abundances in size bin
                    bin_fish <- 1:no_w %in% bin.id
                    bin_fish <- sweep(n, 2, bin_fish, "*")
                    # Calculate number of competitors from all refuge users in bin k
                    refuge_user_idx <- which(params@species_params$refuge_user == TRUE)
                    competitors <- bin_fish[refuge_user_idx, , drop = FALSE] %*% params@dw
                    competitor_density[idx] <- sum(competitors)
                    # Set vulnerability for fish in size bin for this species
                    refuge[i, bin.id] <- ifelse(competitor_density[idx] == 0,
                        max_protect,
                        tau * new_rd[k] / competitor_density[idx]
                    )
                }
            }
        } else if (use_dummy_fish_bins) {
            # Handle non-species-specific bins
            competitor_density <- numeric(length(new_rd))
            for (k in seq_len(nrow(new_rd))) {
                bin.id <- bin_id_list[[k]]
                bin_fish <- 1:no_w %in% bin.id
                bin_fish <- sweep(n, 2, bin_fish, "*")
                competitors <- bin_fish %*% params@dw
                # Remove species that don't use refuge
                sp <- params@species_params$species
                sp <- sp[params@species_params$refuge_user == TRUE]
                competitors <- competitors[sp, ]
                competitor_density[k] <- sum(competitors)
                refuge[, bin.id] <- ifelse(competitor_density[k] == 0,
                    max_protect,
                    tau * new_rd[k] / competitor_density[k]
                )
            }
        }
        # Make sure none of the values are higher than max_protect
        refuge[refuge > max_protect] <- max_protect
        # Account for vulnerability of species that don't utilize refuge
        vulnerable <- 1 - (refuge_user * refuge)
    }

    return(vulnerable)
}


# Old reef vulnerable
reefVulnerable <- function(
  params, n, n_pp, n_other, t,
  new_rd = reefDegrade(params, n, n_pp, n_other, t, ...)
) {
    # Extract relevant data from params
    refuge_params <- params@other_params[["refuge_params"]]
    method_params <- params@other_params[["method_params"]]

    method_params$refuge_density <- new_rd

    # Set parameters used with all methods
    w_settle <- refuge_params$w_settle
    max_protect <- refuge_params$max_protect
    tau <- refuge_params$tau

    # Pull no of species and size bins
    no_w <- length(params@w)
    no_sp <- dim(params@interaction)[1]

    # Store which functional groups use refuge
    refuge_user <- params@species_params$refuge_user

    # Static methods -----------------------------------------------------------
    static <- c("sigmoidal", "binned", "noncomplex")

    if (is.element(refuge_params$method, static)) {
        refuge <- params@other_params$refuge
        vulnerable <- 1 - refuge

        # Competitive method -------------------------------------------------------
    } else if (refuge_params$method == "competitive") {
        # Initialize empty list to hold number of competitors for each bin
        competitor_density <- numeric(length(method_params$refuge_density))

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
            sp <- sp[params@species_params$refuge_user == TRUE]
            competitors <- competitors[sp, ]

            # sum competitors from all species groups for refuge bin k
            competitor_density[k] <- sum(competitors)

            # Set vulnerability for fish in size bin based on the number of
            # available refuges and the number of competitors
            refuge[, bin.id] <- ifelse(competitor_density[k] == 0,
                max_protect,
                tau * method_params$refuge_density[k] / competitor_density[k]
            )
        }
        # Make sure none of the values are higher than max_protect
        refuge[refuge > max_protect] <- max_protect
        # Account for vulnerability of species that don't utilize refuge
        vulnerable <- 1 - (refuge_user * refuge)
    }

    return(vulnerable)
}

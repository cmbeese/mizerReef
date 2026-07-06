# New getRefuge
getRefuge <- function(params, use_dummy_fish_bins = TRUE, ...) {
    # object - Check if mizerParams is valid
    assert_that(is(params, "MizerParams"))

    # Extract relevant data from params
    refuge_params <- params@other_params$refuge_params
    method_params <- params@other_params$refuge_params[["method_params"]]

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
        params@other_params$refuge_params$refuge <- refuge

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
        params@other_params$refuge_params$refuge <- refuge
        params@other_params$refuge_params$refuge_lengths <- refuge_lengths
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
        params@other_params$refuge_params$refuge <- refuge
        params@other_params$refuge_params$bin.id <- bin.id
        # store length bins by functional group in params object
        start_l.i <- t(do.call(rbind, start_l.i))
        end_l.i <- t(do.call(rbind, end_l.i))
        refuge_lengths <- cbind(start_l.i, end_l.i)
        row.names(refuge_lengths) <- sp$species
        params@other_params$refuge_params$refuge_lengths <- refuge_lengths
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
        params@other_params$refuge_params$bin.id <- bin.id
        # Store length bins by species group in a data frame
        start_l.i <- t(do.call(rbind, start_l.i))
        end_l.i <- t(do.call(rbind, end_l.i))
        refuge_lengths <- cbind(start_l.i, end_l.i)
        row.names(refuge_lengths) <- sp$species
        params@other_params$refuge_params$refuge_lengths <- refuge_lengths
        params@time_modified <- lubridate::now()
    }
    return(params)
}


# Old get refuge
getRefuge <- function(params, ...) {
    # object - Check if mizerParams is valid
    assert_that(is(params, "MizerParams"))

    # Extract relevant data from params
    refuge_params <- params@other_params[["refuge_params"]]
    method_params <- params@other_params[["method_params"]]

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
    tau <- refuge_params$tau

    # Store which functional groups use refuge
    refuge_user <- sp$refuge_user

    # Noncomplex reef ----------------------------------------------------------
    if (refuge_params$method == "noncomplex") {
        # Create matrix to store proportions for each species
        refuge <- matrix(0, nrow = no_sp, ncol = no_w)
        rownames(refuge) <- rownames(params@initial_n)
        colnames(refuge) <- colnames(params@initial_n)

        # store refuge and bin indices in params object
        params@other_params[["refuge"]] <- refuge

        # Sigmoidal method ---------------------------------------------------------
    } else if (refuge_params$method == "sigmoidal") {
        # Pull slope and proportion of fish to be protected from method_params
        prop_protect <- method_params$prop_protect
        slope <- method_params$slope

        # Convert length to weight to determine refuge capacity
        W_refuge <- a_bar * method_params$L_refuge^b_bar

        # Calculate sigmoid using threshold weights - no organisms smaller
        # than w_settle or larger than W_refuge can utilize refuge
        denom <- 1 + exp(slope * (w - W_refuge))
        ref <- ifelse(w > w_settle, prop_protect / denom, 0)

        # Make sure none of the values are higher than maximum protection allowed
        ref[ref > max_protect] <- max_protect

        # Create matrix to store proportions for each species
        refuge <- matrix(rep(ref), nrow = no_sp, ncol = no_w, byrow = TRUE)
        rownames(refuge) <- rownames(params@initial_n)
        colnames(refuge) <- colnames(params@initial_n)

        # Account for species that don't utilize refuge
        refuge <- refuge_user * refuge

        # store refuge in params object
        params@other_params[["refuge"]] <- refuge

        # Find L_refuge by species & store in data frame
        L_refuge.i <- (W_refuge / sp[["a"]])^(1 / sp[["b"]])
        refuge_lengths <- data.frame(sp$species, L_refuge.i)
        params@other_params[["refuge_lengths"]] <- refuge_lengths

        # Save time parameters were modified
        params@time_modified <- lubridate::now()

        # Binned method ------------------------------------------------------------
    } else if (refuge_params$method == "binned") {
        # Initialize storage
        ref <- rep(0, no_w)
        start_l.i <- list(1)
        end_l.i <- list(1)
        bin.id <- list(1)
        no_bins <- nrow(method_params)

        # Loop through each refuge bin
        for (k in 1:no_bins) {
            # Calculate start and end of weight bins for a dummy fish
            start_w <- a_bar * method_params$start_L[[k]]^b_bar
            end_w <- a_bar * method_params$end_L[[k]]^b_bar

            # Set threshold weight - no organisms smaller than w_settle
            start_w[start_w < w_settle] <- w_settle

            # Gives indices of fish in size range to protect
            bin.id[[k]] <- which(w >= start_w & w <= end_w)

            # Refuge
            ref[bin.id[[k]]] <- method_params$prop_protect[k]

            # Calculate length bins for each species
            start_l.i[[k]] <- (start_w / sp[["a"]])^(1 / sp[["b"]])
            names(start_l.i)[[k]] <- c(paste("start", k, sep = ""))
            end_l.i[[k]] <- (end_w / sp[["a"]])^(1 / sp[["b"]])
            names(end_l.i)[[k]] <- c(paste("end", k, sep = ""))
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
        params@other_params[["refuge"]] <- refuge
        params@other_params[["bin.id"]] <- bin.id

        # store length bins by functional group in params object
        start_l.i <- t(do.call(rbind, start_l.i))
        end_l.i <- t(do.call(rbind, end_l.i))
        refuge_lengths <- cbind(start_l.i, end_l.i)
        row.names(refuge_lengths) <- sp$species
        params@other_params[["refuge_lengths"]] <- refuge_lengths

        # Save time parameters were modified
        params@time_modified <- lubridate::now()

        ## Competitive method ------------------------------------------------------
    } else if (refuge_params$method == "competitive") {
        # Initialize empty list to hold number of competitors for each bin
        competitor_density <- numeric(nrow(method_params))

        # Empty list to hold indices of fish protected by each bin
        bin.id <- list(1)
        start_l.i <- list(1)
        end_l.i <- list(1)

        # Loop through each refuge bin
        for (k in 1:nrow(method_params)) {
            # Calculate start and end of weight bin k
            start_w <- a_bar * method_params$start_L[[k]]^b_bar
            end_w <- a_bar * method_params$end_L[[k]]^b_bar

            # No organisms smaller than w_settle can use refuge
            start_w[start_w < w_settle] <- w_settle

            # Find indices of fish within size bin k
            bin.id[[k]] <- which(params@w >= start_w & params@w <= end_w)

            # Calculate length bins for each species
            start_l.i[[k]] <- (start_w / sp[["a"]])^(1 / sp[["b"]])
            names(start_l.i)[[k]] <- c(paste("start", k, sep = ""))
            end_l.i[[k]] <- (end_w / sp[["a"]])^(1 / sp[["b"]])
            names(end_l.i)[[k]] <- c(paste("end", k, sep = ""))
        }

        # Store indices of each bin
        params@other_params[["bin.id"]] <- bin.id

        # Store length bins by species group in a data frame
        start_l.i <- t(do.call(rbind, start_l.i))
        end_l.i <- t(do.call(rbind, end_l.i))
        refuge_lengths <- cbind(start_l.i, end_l.i)
        row.names(refuge_lengths) <- sp$species
        params@other_params[["refuge_lengths"]] <- refuge_lengths

        # Save time parameters were modified
        params@time_modified <- lubridate::now()
    }
    return(params)
}

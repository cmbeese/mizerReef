# Load mizer
library(mizer)
library(assertthat)
# Call NS params
params <- NS_params
# Set default a and b values for each species
params <- set_species_param_default(params, 'a', 0.025)
params <- set_species_param_default(params, 'b', 3)
no_sp <- dim(params@interaction)[1]

# Creat some tester complexity data
sig <- data.frame(L_refuge = 15, prop_protect = 0.2)
bin <- data.frame(start_L = seq(0, 45, 5),
                  end_L = seq(5, 50, 5),
                  prop_protect = c(1.0, 0.8, 0.6, 0.4, 0.4,
                                   0.4, 0.3, 0.3,   0,   0))
comp <- data.frame(start_L = seq(0, 45, 5),
                   end_L = seq(5, 50, 5),
                   refuge_density = c(0, 10^5, 10^5, 0, 0,
                                      0, 10^6, 10^3, 0, 0))

# Use set refuge function 
# sig_params <- setRefuge(params, method = 'sigmoidal',
#                         method_params = sig,
#                         refuge_user = rep(TRUE, no_sp),
#                         bad_pred = rep(TRUE, no_sp))

bin_params <- setRefuge(params, method = 'binned',
                        method_params = bin,
                        refuge_user = rep(TRUE, no_sp),
                        bad_pred = rep(TRUE, no_sp))

com_params <- setRefuge(params, method = 'competitive',
                        method_params = bin,
                        refuge_user = rep(TRUE, no_sp),
                        bad_pred = rep(TRUE, no_sp))

params <- bin_params

getRefuge <- function(params, n, n_pp, n_other, t = 0, ...) {
    
    # Extract relevant data from params
    refuge_params <- params@other_params[['refuge_params']]
    method_params <- params@other_params[['method_params']]
    
    # Pull values from params
    w <- params@w
    sp <- params@species_params
    no_w <- length(params@w)
    no_sp <- dim(params@interaction)[1]
    
    # Initialize storage for the array of refuge proportions (profile)
    refuge <- matrix(0, nrow = no_sp, ncol = no_w)
    rownames(refuge) <- rownames(params@initial_n)
    colnames(refuge) <- colnames(params@initial_n)
    
    # Set parameters used with all methods
    a_bar       <- refuge_params$a_bar
    b_bar       <- refuge_params$b_bar
    w_settle    <- refuge_params$w_settle
    max_protect <- refuge_params$max_protect
    tau         <- refuge_params$tau

    # Store which functional groups use refuge
    refuge_user <- sp$refuge_user
    
    # Sigmoidal method ---------------------------------------------------------
    if (refuge_params$method == "sigmoidal"){
        
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
        refuge[refuge > max_protect] = max_protect
        
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
            bin.id[[k]] <- which(params@w >= start_w & params@w <= end_w)
            
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
    
    # Data method --------------------------------------------------------------
    } else if (refuge_params$method == "competitive") {
        
        # Initialize empty list to hold number of competitors for each bin
        competitor_density = numeric(nrow(method_params))
        # Empty list to hold indices of fish protected by each bin
        bin.id = list(1)
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
    


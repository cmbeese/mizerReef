reefVulnerable <- function(params, n, n_pp, n_other, t = 0, ...) {
    
    # Extract relevant data from params
    refuge_params <- params@other_params[['refuge_params']]
    method_params <- params@other_params[['method_params']]
    
    # Set parameters used with all methods
    w_settle    <- refuge_params$w_settle
    max_protect <- refuge_params$max_protect
    tau         <- refuge_params$tau
    
    # Pull no of spcies and size bins
    no_w <- length(params@w)
    no_sp <- dim(params@interaction)[1]
    
    # Store which functional groups use refuge
    refuge_user <- params@species_params$refuge_user
    
    # Static methods -----------------------------------------------------------
    static = c("sigmoidal", "binned")
    
    if (is.element(refuge_params$method, static)){
        
        refuge <- params@other_params$refuge
        vulnerable <- 1 - refuge
    
    # Competitive method -------------------------------------------------------
    } else if (refuge_params$method == "competitive") {
        
        # Initialize empty list to hold number of competitors for each bin
        competitor_density = numeric(nrow(method_params))
        
        # Initialize storage for the array of refuge proportions
        refuge <- matrix(0, nrow = no_sp, ncol = no_w)
        rownames(refuge) <- rownames(params@initial_n)
        colnames(refuge) <- colnames(params@initial_n)
        
        # Loop through each refuge bin
        for (k in 1:nrow(method_params)) {
            
            k = 2
            # Get indices of fish in size bin k
            bin.id <- params@other_params$bin.id[[k]]
            
            # Creat vector of zeroes and ones
            bin_fish <- 1:no_w %in% bin.id
            
            # Calculate number of competitors from each functional group in bin k
            competitors <- (n * bin_fish) %*% params@dw
            
            # Eliminate functional groups that don't use refuge and sum
            competitor_density[k] <- sum(refuge_user * competitors)
            
            # Set vulnerability for fish in size bin based on the number of
            # available refuges and the number of competitors
            refuge[,bin.id] <- ifelse(competitor_density[k] == 0, 
                                            max_protect,
                    tau * method_params$refuge_density[k]/competitor_density[k])
            
            # Make sure none of the values are higher than maximum protection allowed
            refuge[refuge > max_protect] = max_protect
            
            # Account for species that don't utilize refuge
            vulnerable = 1 - (refuge_user*refuge)
        }
    }
    return(vulnerable)
}

# Testing new function ---------------------------------------------------------
library(mizerReef)
library(mizer)

    # params <- bonaire_model
    params <- karpata_model
    # params <- NS_params
    vulnerable <- getVulnerable(params)
    pred_rate <- getPredRate(params)

    # Function steps
        no_sp <- nrow(params@species_params)
        no_w <- length(params@w)
        no_w_full <- length(params@w_full)
        
        # Find number of size bins in resource spectra smaller than smallest fish
        p <- no_w_full - no_w  
        
        # Add columns for entire model size range to vulnerability matrix
        resource_vul <- matrix(1, nrow = no_sp, ncol = p)
        
        # Add resource vulnerability to vulnerability matrix
        vul_pp <- cbind(resource_vul, vulnerable)
        
        # Find indices of predator species whose foraging is hindered by refuge
        bad_pred  <- which(params@species_params$bad_pred == TRUE)
        good_pred <- which(params@species_params$bad_pred == FALSE)
        
        # Create list of vulnerabilities for each predator
        vul <- vector("list", no_sp)
        vul[bad_pred] <- list(vul_pp)
        vul[good_pred] <- list(matrix(1, nrow = no_sp, ncol = ncol(vul_pp)))
        
        # Loop through predator species to calculate predation mortality on
        # each prey species & size by all predators
        pm <- matrix(0, no_sp, length(params@w_full))
        
        for (i in 1:no_sp){
            i = 1
            # Vulnerability rate of all prey, including resource
            # (species by size) to predator i
            v <- vul[[i]]
            # Predation rate of predator species i on all prey (by size)
            pr_i <- pred_rate[i,]
            pr_i_mat <- matrix(rep(pr_i, each = nrow(v)), 
                           nrow = nrow(v), ncol = ncol(v))
            # vul*pr_i predation mortality on prey (species by size) by 
            # predator i, sum across all predators by adding pm
            pm <- pm + v*pr_i_mat
        }

        # Get index of species that have grown out of the resource spectrum
        idx_sp <- (length(params@w_full) - 
                       length(params@w) + 1):length(params@w_full)
        
        # Account for interaction of species
        pred_mort <- base::t(params@interaction) %*% pm[, idx_sp, drop = FALSE]
        
        # Plot to check
        x <- log10(params@w)
        
        plot(x, pred_mort[1,])
        plot(x, pred_mort[2,])
        plot(x, pred_mort[3,])
        plot(x, pred_mort[5,])
        plot(x, pred_mort[6,])
        plot(x, pred_mort[7,])
        plot(x, pred_mort[8,])
        plot(x, pred_mort[9,])
        plot(x, pred_mort[10,])
        
        
        
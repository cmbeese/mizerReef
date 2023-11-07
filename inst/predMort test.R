# 

# test vulnerability
vulnerable = matrix(data = 0.5, nrow = 5, ncol = 10)

# test pred rate
pred_rate = matrix(data = rep(1:5), nrow = 5, ncol = 10)


# Find indices of fish that have grown out of the resource spectrum
idx_sp <- (length(params@w_full) - 
               length(params@w) + 1):length(params@w_full)
no_sp <- nrow(params@species_params)

# Find indices of predator species whose foraging is hindered by refuge
bad_pred  <- which(params@species_params$bad_pred == TRUE)
good_pred <- which(params@species_params$bad_pred == FALSE)

# Create list of vulnerabilities for each predator
vul <- vector("list", no_sp)
vul[bad_pred] <- list(vulnerable)
vul[good_pred] <- list(1)
p <- matrix(0, dim(n), dim(n))
i = 2

p[i,] <- vul[[i]] * pred_rate[i, idx_sp, drop = FALSE]

v <- vul[[i]]
pr_i <- pred_rate[i, idx_sp, drop = TRUE]

# Loop through predator species
p <- matrix(0, dim(n), dim(n))
for (i in 1:no_sp){
    vul <- vu
    pm_by_prey <- lapply(vul[[i]], function(row) row * row_to_multiply)
    p <-  * pred_rate[, idx_sp, drop = FALSE]
}

return(base::t(params@interaction) %*% p)
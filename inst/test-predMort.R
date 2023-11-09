
# test vulnerability
vulnerable <- matrix(data = 0.5, nrow = 5, ncol = 10)

# test pred rate
pred_rate <- matrix(data = rep(1:5), nrow = 5, ncol = 15)

# Find indices of fish that have grown out of the resource spectrum
idx_sp <- 6:15
no_sp <- 5

# Find indices of predator species whose foraging is hindered by refuge
bad_pred  <- c(1,3,5)
good_pred <- c(2,4)

# Two dimensional method with loop
pm1 <- function(vulnerable, 
                pred_rate, 
                bad_pred, 
                good_pred, 
                idx_sp,
                no_sp){
    # Create list of vulnerabilities for each predator
    vul <- vector("list", no_sp)
    vul[bad_pred] <- list(vulnerable)
    vul[good_pred] <- list(matrix(1, nrow = no_sp, ncol = length(idx_sp)))
    
    # Loop through predator species
    pm <- matrix(0, no_sp, length(idx_sp))
    for (i in 1:no_sp){
        v <- vul[[i]]
        pr_i <- pred_rate[i, idx_sp]
        pm <- pm + v*pr_i
    }
}

# Two dimensional method with lapply
pm2 <- function(vulnerable, 
                pred_rate, 
                bad_pred, 
                good_pred, 
                idx_sp,
                no_sp){
    
    # Create list of vulnerabilities for each predator
    vul <- vector("list", no_sp)
    vul[bad_pred] <- list(vulnerable)
    vul[good_pred] <- list(matrix(1, nrow = no_sp, ncol = length(idx_sp)))

    # Calculate the product of vulnerabilities and predator rates using lapply
    pm2_list <- lapply(1:no_sp, function(i) vul[[i]] * pred_rate[i, idx_sp])
    pm2 <- Reduce('+', pm2_list)
    
}

# Three dimensional method with lapply
pm3 <- function(vulnerable, 
                pred_rate, 
                bad_pred, 
                good_pred, 
                idx_sp,
                no_sp){
    
    # Create a three-dimensional matrix for vulnerabilities
    vul <- array(1, dim = c(no_sp, ncol(vulnerable), no_sp))
    vul[,,bad_pred] <- vulnerable
    vul[,,good_pred] <- matrix(1, nrow = no_sp, ncol = length(idx_sp))
    
    pm3_list <- lapply(1:no_sp, function(i) vul[,,i] * pred_rate[i, idx_sp])
    pm3 <- Reduce('+', pm3_list)
}

# Speed test 
library(microbenchmark)
results <- microbenchmark(
    pm1(vulnerable, pred_rate, bad_pred, good_pred, idx_sp, no_sp),
    pm2(vulnerable, pred_rate, bad_pred, good_pred, idx_sp, no_sp),
    pm3(vulnerable, pred_rate, bad_pred, good_pred, idx_sp, no_sp),
    times = 10000 
)

# Pm1 - loops method is winner!

#return(base::t(params@interaction) %*% p)
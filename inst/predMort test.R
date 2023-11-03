# Load species parameter data
cbn_species <- read.csv("data/cbn_species.csv", row.names = 1, header = TRUE)
cbn_int     <- read.csv("data/cbn_int.csv",     row.names = 1, header = TRUE)
cbn_UR_int  <- read.csv("data/cbn_UR_int.csv",  row.names = 1, header = TRUE)

# Create some tester refuge scenarios
methods <- c("simple", "binned", "data")
beese_2023  <- data.frame(start_L = seq(0, 45, 5),
                          end_L = seq(5, 50, 5),
                          prop_protect = c(1, 1, 0.9, rep(0,7)))
## SET MODEL -------------------------------------------------------------------
params <- newReefParams(species_params = cbn_species,
                        interaction = cbn_int,
                        UR_interaction = cbn_UR_int,
                        method = methods[2],
                        method_params = beese_2023)
# Get vulnerability
vulnerable = getVulnerable(params)
pred_rate = getPredRate(params)
n = params@initial_n

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
i = 1

p <- vul[[i]] * pred_rate[, idx_sp, drop = FALSE]

# Loop through predator species
p <- matrix(0, dim(n), dim(n))
for (i in 1:no_sp){
    p <- vul[[i]] * pred_rate[, idx_sp, drop = FALSE]
}

return(base::t(params@interaction) %*% p)
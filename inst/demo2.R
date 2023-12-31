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
                        method_params = comp,
                        refuge_user = rep(TRUE, no_sp),
                        bad_pred = rep(TRUE, no_sp))
library(here)

# Rubble trajectory
rubble_scale <- read.csv(here("inst/rubble.csv"), header = FALSE)
rubble_scale <- as.matrix(rubble_scale)

save(rubble_scale,   file = "data/rubble_scale.rda")

# Algae trajectory
algae_scale <- read.csv(here("inst/algae.csv"), header = FALSE)
algae_scale <- as.matrix(algae_scale)

save(algae_scale,   file = "data/algae_scale.rda")

# Recovery trajectory
recovery_scale <- read.csv(here("inst/recovery.csv"), header = FALSE)
recovery_scale <- as.matrix(recovery_scale)

save(recovery_scale,   file = "data/recovery_scale.rda")

# No degradation - for testing
constant_scale <- matrix(1, nrow = nrow(algae_scale), ncol = ncol(algae_scale))
save(constant_scale,   file = "data/constant_scale.rda")

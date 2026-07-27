# Regenerate the built-in degradation trajectory scaling matrices
# (rubble_scale, algae_scale, recovery_scale) from inst/data-csv/deg_scales.csv
#
# Each trajectory in deg_scales.csv stores a *change* (delta) in refuge
# density per refuge size bin (rows) and year (columns 1-15); the scaling
# matrices used by setDegradation()/reefDegrade() are multipliers, i.e.
# 1 + delta, applied to the previous year's refuge density.
library(here)

deg_scales <- read.csv(here("inst/data-csv/deg_scales.csv"))
year_cols <- as.character(1:15)

make_scale_matrix <- function(trajectory) {
    d <- deg_scales[deg_scales$trajectory == trajectory, ]
    m <- 1 + as.matrix(d[, paste0("X", year_cols)])
    dimnames(m) <- list(d$ref_size, year_cols)
    m
}

# Rubble trajectory
rubble_scale <- make_scale_matrix("rubble")
save(rubble_scale, file = "data/rubble_scale.rda")

# Algae trajectory
algae_scale <- make_scale_matrix("algae")
save(algae_scale, file = "data/algae_scale.rda")

# Recovery trajectory
recovery_scale <- make_scale_matrix("recovery")
save(recovery_scale, file = "data/recovery_scale.rda")

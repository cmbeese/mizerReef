
## Setup - load packages ----------------------------------
library(ggplot2)
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(assertthat)
library(here)

params <- karpata_model


# Normal
sim1 <- project(params, tmax = 10) # works fine
b <- plotBiomass(sim1, start_time = 1, end_time = 10, total = TRUE, 
                 return_data = TRUE)
p <- plotProductivity(sim1, start_time = 1, end_time = 10, 
                      return_data = TRUE) # also works
df <- mizer::plotBiomass(sim1,total = TRUE, return_data = TRUE)

# Add carrying capacity
params2 <- setURcapacity(params)
params2 <- reefSteady(params2)
sim2 <- project(params2, tmax = 10) # this also works
b <- plotBiomass(sim2, start_time = 1, end_time = 10)
p <- plotProductivity(sim2, start_time = 1, end_time = 10)

# Add degradation without carrying capacity
params3 <- setDegradation(params, trajectory = "rubble", 
                          deg_scale = constant_scale)
params3 <- reefSteady(params3)
sim3 <- project(params3, tmax = 10)
b <- plotBiomass(sim3, start_time = 1, end_time = 10)
p <- plotProductivity(sim3, start_time = 1, end_time = 10)
print(b)
print(p)

# Add degradation with carrying capacity
params4 <- setDegradation(params, trajectory = "rubble", 
                          deg_scale = rubble_scale)
params4 <- setURcapacity(params4)
params4 <- reefSteady(params4)
sim4 <- project(params4, tmax = 10)
b <- plotBiomass(sim4, start_time = 1, end_time = 10)
p <- plotProductivity(sim4, start_time = 1, end_time = 10)
print(b)
print(p)

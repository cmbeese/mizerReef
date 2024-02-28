
## Setup - load packages ----------------------------------
library(ggplot2)
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(assertthat)
library(here)

params <- karpata_model


# Normal
sim1 <- project(params, t_max = 30) # works fine
b <- plotBiomass(sim1, total = TRUE,  return_data = TRUE)
p <- plotProductivity(sim1, return_data = TRUE) # also works
print(plotBiomass(sim1))

# Add carrying capacity
params2 <- setURcapacity(params)
params2 <- reefSteady(params2)
sim2 <- project(params2, t_max = 30) # this also works
b <- plotBiomass(sim2)
p <- plotProductivity(sim2)
print(b)
print(p)

# Add degradation without carrying capacity
params3 <- setDegradation(params, trajectory = "rubble", 
                          deg_scale = constant_scale)
params3 <- reefSteady(params3)
sim3 <- project(params3, t_max = 30)
b <- plotBiomass(sim3)
p <- plotProductivity(sim3)
print(b)
print(p)

# Add degradation with carrying capacity
params4 <- setDegradation(params, trajectory = "rubble", 
                          deg_scale = rubble_scale)
params4 <- setURcapacity(params4)
params4 <- reefSteady(params4)
params4 <- reefSteady(params4)
sim4 <- project(params4, t_max = 30)
b <- plotBiomass(sim4)
p <- plotProductivity(sim4)
print(b)
print(p)

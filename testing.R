library(tidyverse)
library(ggpubr)

# Reading the data
projections <- read.csv("anomaly_df.csv")
SST_dev <- projections[['anomaly']]

# Initializing parameters
patch_area <- c(1)
number_patches <- length(patch_area)
timesteps <- length(SST_dev)
r <- 0.3
K <- 101.3

# Initializing arrays
population <- array(NA, dim = c(timesteps, number_patches))
population[1, ] <- c(10)
harvest <- array(NA, dim = c(timesteps, number_patches))

patch_area_sequences <- seq(1, 1, by = 1)
fishing_effort_sequences <- seq(0, 1, by = 0.1)

parameter_grid <- expand.grid(patch_area = patch_area_sequences,
                              fishing_effort = fishing_effort_sequences)

# Saving outputs
outcome_population <- matrix(0, ncol = 1, nrow = nrow(parameter_grid))
outcome_harvest <- matrix(0, ncol = 1, nrow = nrow(parameter_grid))

# Simulation loop
for (iter in 1:nrow(parameter_grid)){
  for (t in 2:timesteps) {
    population[t, ] <- population[t-1, ] + r * population[t-1, ] * (1 - population[t-1, ] / K)
    harvest[t, ] <- population[t, ] * (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]])))
    population[t, ] <- population[t, ] * (1 - (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]))))
  }
  outcome_population[iter] <- population[t, ] #mean(population[(t-19):t])
  outcome_harvest[iter] <- sum(harvest[, ]) # sum of harvest across all patches
}

# Combining results with parameter grid
outcome_harvest <- cbind(parameter_grid, outcome_harvest)

# Plotting
ggplot(outcome_harvest, aes(x = fishing_effort, y = outcome_harvest)) +
  geom_line() +
  labs(x = "Fishing Effort", y = "Total Harvest") +
  theme_minimal()










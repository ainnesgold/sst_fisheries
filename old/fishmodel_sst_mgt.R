source('dispersal.R')
S <- 0.5
patch_area <- c(0.8, 0.2)
number_patches <- length(patch_area)

tmp <- dispersal(patch_area, number_patches, S)
#tmp1 <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, S)

fishing_effort = 0.1
catchability = 1

timesteps = 80
r = 0.3
K = 101.3

projections <- read.csv("projections.csv")
SST_dev <- projections[['anomaly']]
SST_dev_higheropttemp <- projections[['anomaly_higheropttemp']]
SST_dev_loweropttemp <- projections[['anomaly_loweropttemp']]

mean_temps_crop <- read.csv("mean_temps_crop.csv")
mean_temps_series <- mean_temps_crop[['extracted_mean_temps']]

population <- array(NA, dim = c(timesteps, number_patches))
population[1,] <- c(10,10)
harvest <- array(NA, dim = c(timesteps, number_patches))

population_r1 <- array(NA, dim = c(timesteps, number_patches))
population_r1[1,] <- c(10,10)
harvest_r1 <- array(NA, dim = c(timesteps, number_patches))

population_r2 <- array(NA, dim = c(timesteps, number_patches))
population_r2[1,] <- c(10,10)
harvest_r2 <- array(NA, dim = c(timesteps, number_patches))

population_r3 <- array(NA, dim = c(timesteps, number_patches))
population_r3[1,] <- c(10,10)
harvest_r3 <- array(NA, dim = c(timesteps, number_patches))

population_K1 <- array(NA, dim = c(timesteps, number_patches))
population_K1[1,] <- c(10,10)
harvest_K1 <- array(NA, dim = c(timesteps, number_patches))

population_K2 <- array(NA, dim = c(timesteps, number_patches))
population_K2[1,] <- c(10,10)
harvest_K2 <- array(NA, dim = c(timesteps, number_patches))


r_temp_1 <- array(NA, dim = c(timesteps, 1))
r_temp_2 <- array(NA, dim = c(timesteps, 1))
r_temp_3 <- array(NA, dim = c(timesteps, 1))

K_temp_1 <- array(NA, dim = c(timesteps, 1))
K_temp_2 <- array(NA, dim = c(timesteps, 1))




for (t in 2:timesteps) {
  
    #no temp
    population[t,] <- population[t-1,] + r * population[t-1,] * (1 - population[t-1,] / K)
    harvest[t,] <- population[t,] * (1 - exp(-fishing_effort * catchability)) * patch_area
    population[t,] <- population[t,] * (1 - (1 - exp(-fishing_effort * catchability)))
    population[t,] <- population[t,] %*% tmp
    
    #temp dependent r - "current" temp as optimal
    r_temp_1[t] <- 0.3 + 0*SST_dev[t] + -0.0037*SST_dev[t]^2
    population_r1[t,] <- population_r1[t-1,] + r_temp_1[t] * population_r1[t-1,] * (1 - population_r1[t-1,] / K)
    harvest_r1[t,] <- population_r1[t,] * (1 - exp(-fishing_effort * catchability)) * patch_area
    population_r1[t,] <- population_r1[t,] * (1 - (1 - exp(-fishing_effort * catchability)))
    population_r1[t,] <- population_r1[t,] %*% tmp
    
    #temp dependent r - higher optimal temp
    r_temp_2[t] <- a + b*SST_dev_higheropttemp[t] + c*SST_dev_higheropttemp[t]^2
    population_r2[t,] <- population_r2[t-1,] + r_temp_2[t] * population_r2[t-1,] * (1 - population_r2[t-1,] / K)
    harvest_r2[t,] <- population_r2[t,] * (1 - exp(-fishing_effort * catchability)) * patch_area
    population_r2[t,] <- population_r2[t,] * (1 - (1 - exp(-fishing_effort * catchability)))
    population_r2[t,] <- population_r2[t,] %*% tmp
    
    #temp dependent r - lower optimal temp
    r_temp_3[t] <- a + b*SST_dev_loweropttemp[t] + c*SST_dev_loweropttemp[t]^2
    population_r3[t,] <- population_r3[t-1,] + r_temp_3[t] * population_r3[t-1,] * (1 - population_r3[t-1,] / K)
    harvest_r3[t,] <- population_r3[t,] * (1 - exp(-fishing_effort * catchability)) * patch_area
    population_r3[t,] <- population_r3[t,] * (1 - (1 - exp(-fishing_effort * catchability)))
    population_r3[t,] <- population_r3[t,] %*% tmp
    
    #temp dependent K - piece wise linear function
    K_temp_1[t] <- -2.1408 * mean_temps_series[t] + 160.253358
    if (K_temp_1[t] < 10) {
      K_temp_1[t] <- 10
    } else if (K_temp_1[t] > 101.3) {
      K_temp_1[t] <- 101.3
    }
    population_K1[t,] <- population_K1[t-1,] + r * population_K1[t-1,] * (1 - population_K1[t-1,] / K_temp_1[t])
    harvest_K1[t,] <- population_K1[t,] * (1 - exp(-fishing_effort * catchability)) * patch_area
    population_K1[t,] <- population_K1[t,] * (1 - (1 - exp(-fishing_effort * catchability)))
    population_K1[t,] <- population_K1[t,] %*% tmp
    
    #temp dependent K - quadratic function
    K_temp_2[t] <- 101.3 + 0*SST_dev[t] + -0.7*SST_dev[t]^2
    if (K_temp_2[t] < 10) {
      K_temp_2[t] <- 10 }
    else if (K_temp_2[t] > 101.3) {
      K_temp_2[t] <- 101.3
    }
    population_K2[t,] <- population_K2[t-1,] + r * population_K2[t-1,] * (1 - population_K2[t-1,] / K_temp_2[t])
    harvest_K2[t,] <- population_K2[t,] * (1 - exp(-fishing_effort * catchability)) * patch_area
    population_K2[t,] <- population_K2[t,] * (1 - (1 - exp(-fishing_effort * catchability)))
    population_K2[t,] <- population_K2[t,] %*% tmp
  
}

plot(population[,1])

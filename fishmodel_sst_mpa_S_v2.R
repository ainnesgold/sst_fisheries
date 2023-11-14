library(tidyverse)
library(ggpubr)
#Varying S, sst, and MPA
#not looking very interesting at the moment

source('dispersal.R')

projections <- read.csv("anomaly_df.csv")
SST_dev <- projections[['anomaly']]
SST_dev_loweropttemp <- projections[['anomaly']] + 1 #larger anomaly = lower optimal temp
SST_dev_higheropttemp <- projections[['anomaly']] - 1 #smaller anomaly = higher optimal temp

timesteps <- length(SST_dev)
r = 0.3
K = c(101.3, 101.3)
S <- seq(0, 1, by = 0.1)

patch_area_sequences <- list(seq(0, 1, by = 0.1))
patch_area_grid <- do.call(expand.grid, patch_area_sequences)
patch_area_grid$Var2 <- 1 - patch_area_grid$Var1
patch_area_list <- split(patch_area_grid, 1:nrow(patch_area_grid))
number_patches <- ncol(patch_area_grid)

fishing_effort <- c(0.2, 0)
catchability <- 1

parameter_grid <- expand.grid(patch_area = patch_area_list,
                              S = S)

#saving outputs
outcome_population <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))

outcome_population_r1 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest_r1 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))

outcome_population_r2 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest_r2 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))

outcome_population_r3 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest_r3 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))

outcome_population_K1 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest_K1 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))

outcome_population_K2 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest_K2 <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))


population <- array(NA, dim = c(timesteps, number_patches))
harvest <- array(NA, dim = c(timesteps, number_patches))

population_r1 <- array(NA, dim = c(timesteps, number_patches))
harvest_r1 <- array(NA, dim = c(timesteps, number_patches))

population_r2 <- array(NA, dim = c(timesteps, number_patches))
harvest_r2 <- array(NA, dim = c(timesteps, number_patches))

population_r3 <- array(NA, dim = c(timesteps, number_patches))
harvest_r3 <- array(NA, dim = c(timesteps, number_patches))

population_K1 <- array(NA, dim = c(timesteps, number_patches))
harvest_K1 <- array(NA, dim = c(timesteps, number_patches))

population_K2 <- array(NA, dim = c(timesteps, number_patches))
harvest_K2 <- array(NA, dim = c(timesteps, number_patches))


r_temp_1 <- array(NA, dim = c(timesteps, 1))
r_temp_2 <- array(NA, dim = c(timesteps, 1))
r_temp_3 <- array(NA, dim = c(timesteps, 1))

K_temp_1 <- array(NA, dim = c(timesteps, 1))
K_temp_2 <- array(NA, dim = c(timesteps, 1))



for (iter in 1:nrow(parameter_grid)){
  
  #dispersal matrix
  tmp <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, parameter_grid[iter, "S"])
  
  #starting numbers
  population[1,] <- c(10,10)
  population_r1[1,] <- c(10,10)
  population_r2[1,] <- c(10,10)
  population_r3[1,] <- c(10,10)
  population_K1[1,] <- c(10,10)
  population_K2[1,] <- c(10,10)
  
  for (t in 2:timesteps) {
    
    #no temp
    if(as.numeric(parameter_grid[['patch_area']][[iter]][1]) == 0 | as.numeric(parameter_grid[['patch_area']][[iter]][2]) == 0){
      patch_zero     <- which(as.numeric(parameter_grid[['patch_area']][[iter]]) == 0) # find which patch is above K
      population[t-1, patch_zero]     <- 0 # force patch above K to equal K
    }
    population[t,] <- population[t-1,] + r * population[t-1,] * (1 - population[t-1,] / K)
    harvest[t,] <- population[t,] * (1 - exp(-fishing_effort * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population[t,] <- population[t,] * (1 - (1 - exp(-fishing_effort * catchability)))
    population[t,] <- population[t,] %*% tmp
    
    #if its above K after dispersal, set it back down to K
    if(population[t, 1] > K[1] & population[t, 2] > K[2]){
      population[t, ] <- K
    }
    else if(population[t, 1] > K[1] | population[t, 2] > K[2]){
      patch_above     <- which(population[t, ] > K) # find which patch is above K
      patch_not_above <- which(!(population[t, ] > K)) # find patch not above K
      spillover       <- population[t, patch_above] -  K[patch_above] # set spillover to the difference between population and K
      
      population[t, patch_above]     <- K[patch_above] # force patch above K to equal K
      population[t, patch_not_above] <- population[t, patch_not_above] + spillover # set the other patch equal to population plus spillover
      
      if(population[t, patch_not_above] > K[patch_not_above]){ # need to check and make sure spillover to the other patch does not push the population over carrying capactiy
        population[t, patch_not_above] <- K[patch_not_above] # if it does then set that patch to carrying capacity after spillover
      }
    } 
    
    
    #temp dependent r - "current" temp as optimal
    r_temp_1[t] <- 0.3 + 0*SST_dev[t] + -0.0037*SST_dev[t]^2
    
    if(as.numeric(parameter_grid[['patch_area']][[iter]][1]) == 0 | as.numeric(parameter_grid[['patch_area']][[iter]][2]) == 0){
      patch_zero     <- which(as.numeric(parameter_grid[['patch_area']][[iter]]) == 0) # find which patch is above K
      population_r1[t-1, patch_zero]     <- 0 # force patch above K to equal K
    }
    population_r1[t,] <- population_r1[t-1,] + r_temp_1[t] * population_r1[t-1,] * (1 - population_r1[t-1,] / K)
    harvest_r1[t,] <- population_r1[t,] * (1 - exp(-fishing_effort * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population_r1[t,] <- population_r1[t,] * (1 - (1 - exp(-fishing_effort * catchability)))
    population_r1[t,] <- population_r1[t,] %*% tmp
    
    if(population_r1[t, 1] > K[1] & population_r1[t, 2] > K[2]){
      population_r1[t, ] <- K
    }
    else if(population_r1[t, 1] > K[1] | population_r1[t, 2] > K[2]){
      patch_above     <- which(population_r1[t, ] > K) # find which patch is above K
      patch_not_above <- which(!(population_r1[t, ] > K)) # find patch not above K
      spillover       <- population_r1[t, patch_above] -  K[patch_above] # set spillover to the difference between population and K
      
      population_r1[t, patch_above]     <- K[patch_above] # force patch above K to equal K
      population_r1[t, patch_not_above] <- population_r1[t, patch_not_above] + spillover # set the other patch equal to population plus spillover
      
      if(population_r1[t, patch_not_above] > K[patch_not_above]){ # need to check and make sure spillover to the other patch does not push the population over carrying capactiy
        population_r1[t, patch_not_above] <- K[patch_not_above] # if it does then set that patch to carrying capacity after spillover
      }
    } 
    
    
    #temp dependent r - higher optimal temp
    r_temp_2[t] <- 0.3 + 0*SST_dev_higheropttemp[t] + -0.0037*SST_dev_higheropttemp[t]^2
    
    if(as.numeric(parameter_grid[['patch_area']][[iter]][1]) == 0 | as.numeric(parameter_grid[['patch_area']][[iter]][2]) == 0){
      patch_zero     <- which(as.numeric(parameter_grid[['patch_area']][[iter]]) == 0) # find which patch is above K
      population_r2[t-1, patch_zero]     <- 0 # force patch above K to equal K
    }
    population_r2[t,] <- population_r2[t-1,] + r_temp_2[t] * population_r2[t-1,] * (1 - population_r2[t-1,] / K)
    harvest_r2[t,] <- population_r2[t,] * (1 - exp(-fishing_effort * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population_r2[t,] <- population_r2[t,] * (1 - (1 - exp(-fishing_effort * catchability)))
    population_r2[t,] <- population_r2[t,] %*% tmp
    
    if(population_r2[t, 1] > K[1] & population_r2[t, 2] > K[2]){
      population_r2[t, ] <- K
    }
    else if(population_r2[t, 1] > K[1] | population_r2[t, 2] > K[2]){
      patch_above     <- which(population_r2[t, ] > K) # find which patch is above K
      patch_not_above <- which(!(population_r2[t, ] > K)) # find patch not above K
      spillover       <- population_r2[t, patch_above] -  K[patch_above] # set spillover to the difference between population and K
      
      population_r2[t, patch_above]     <- K[patch_above] # force patch above K to equal K
      population_r2[t, patch_not_above] <- population_r2[t, patch_not_above] + spillover # set the other patch equal to population plus spillover
      
      if(population_r2[t, patch_not_above] > K[patch_not_above]){ # need to check and make sure spillover to the other patch does not push the population over carrying capactiy
        population_r2[t, patch_not_above] <- K[patch_not_above] # if it does then set that patch to carrying capacity after spillover
      }
    } 
    
    
    #temp dependent r - lower optimal temp
    r_temp_3[t] <- 0.3 + 0*SST_dev_loweropttemp[t] + -0.0037*SST_dev_loweropttemp[t]^2
    
    if(as.numeric(parameter_grid[['patch_area']][[iter]][1]) == 0 | as.numeric(parameter_grid[['patch_area']][[iter]][2]) == 0){
      patch_zero     <- which(as.numeric(parameter_grid[['patch_area']][[iter]]) == 0) 
      population_r3[t-1, patch_zero]     <- 0
    }
    population_r3[t,] <- population_r3[t-1,] + r_temp_3[t] * population_r3[t-1,] * (1 - population_r3[t-1,] / K)
    harvest_r3[t,] <- population_r3[t,] * (1 - exp(-fishing_effort * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population_r3[t,] <- population_r3[t,] * (1 - (1 - exp(-fishing_effort * catchability)))
    population_r3[t,] <- population_r3[t,] %*% tmp
    
    if(population_r3[t, 1] > K[1] & population_r3[t, 2] > K[2]){
      population_r3[t, ] <- K
    }
    else if(population_r3[t, 1] > K[1] | population_r3[t, 2] > K[2]){
      patch_above     <- which(population_r3[t, ] > K) # find which patch is above K
      patch_not_above <- which(!(population_r3[t, ] > K)) # find patch not above K
      spillover       <- population_r3[t, patch_above] -  K[patch_above] # set spillover to the difference between population and K
      
      population_r3[t, patch_above]     <- K[patch_above] # force patch above K to equal K
      population_r3[t, patch_not_above] <- population_r3[t, patch_not_above] + spillover # set the other patch equal to population plus spillover
      
      if(population_r3[t, patch_not_above] > K[patch_not_above]){ # need to check and make sure spillover to the other patch does not push the population over carrying capactiy
        population_r3[t, patch_not_above] <- K[patch_not_above] # if it does then set that patch to carrying capacity after spillover
      }
    } 
    
    
    #temp dependent K - piece wise linear function
    K_temp_1[t] <- -4.95243768 * SST_dev[t] + 101.3
    if (K_temp_1[t] < 10) {
      K_temp_1[t] <- 10
    } else if (K_temp_1[t] > 101.3) {
      K_temp_1[t] <- 101.3
    }
    
    if(as.numeric(parameter_grid[['patch_area']][[iter]][1]) == 0 | as.numeric(parameter_grid[['patch_area']][[iter]][2]) == 0){
      patch_zero     <- which(as.numeric(parameter_grid[['patch_area']][[iter]]) == 0) 
      population_K1[t-1, patch_zero]     <- 0
    }
    population_K1[t,] <- population_K1[t-1,] + r * population_K1[t-1,] * (1 - population_K1[t-1,] / K_temp_1[t])
    harvest_K1[t,] <- population_K1[t,] * (1 - exp(-fishing_effort * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population_K1[t,] <- population_K1[t,] * (1 - (1 - exp(-fishing_effort * catchability)))
    population_K1[t,] <- population_K1[t,] %*% tmp
    
    if(population_K1[t, 1] > K_temp_1[t] & population_K1[t, 2] > K_temp_1[t]){
      population_K1[t, ] <- K_temp_1[t]
    }
    else if(population_K1[t, 1] > K_temp_1[t] | population_K1[t, 2] > K_temp_1[t]){
      patch_above     <- which(population_K1[t, ] > K_temp_1[t]) # find which patch is above K
      patch_not_above <- which(!(population_K1[t, ] > K_temp_1[t])) # find patch not above K
      spillover       <- population_K1[t, patch_above] - K_temp_1[t] # set spillover to the difference between population and K
      
      population_K1[t, patch_above]     <- K_temp_1[t] # force patch above K to equal K
      population_K1[t, patch_not_above] <- population_K1[t, patch_not_above] + spillover # set the other patch equal to population plus spillover
      
      if(population_K1[t, patch_not_above] > K_temp_1[t]){ # need to check and make sure spillover to the other patch does not push the population over carrying capactiy
        population_K1[t, patch_not_above] <- K_temp_1[t] # if it does then set that patch to carrying capacity after spillover
      }
    }
    
    
    
    
    #temp dependent K - quadratic function
    K_temp_2[t] <- 101.3 + 0*SST_dev[t] + -0.7*SST_dev[t]^2
    if (K_temp_2[t] < 10) {
      K_temp_2[t] <- 10 }
    else if (K_temp_2[t] > 101.3) {
      K_temp_2[t] <- 101.3
    }
    
    if(as.numeric(parameter_grid[['patch_area']][[iter]][1]) == 0 | as.numeric(parameter_grid[['patch_area']][[iter]][2]) == 0){
      patch_zero     <- which(as.numeric(parameter_grid[['patch_area']][[iter]]) == 0)
      population_K2[t-1, patch_zero]     <- 0 
    }
    population_K2[t,] <- population_K2[t-1,] + r * population_K2[t-1,] * (1 - population_K2[t-1,] / K_temp_2[t])
    harvest_K2[t,] <- population_K2[t,] * (1 - exp(-fishing_effort * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population_K2[t,] <- population_K2[t,] * (1 - (1 - exp(-fishing_effort * catchability)))
    population_K2[t,] <- population_K2[t,] %*% tmp
    
    if(population_K2[t, 1] > K_temp_2[t] & population_K2[t, 2] > K_temp_2[t]){
      population_K2[t, ] <- K_temp_2[t]
    }
    else if(population_K2[t, 1] > K_temp_2[t] | population_K2[t, 2] > K_temp_2[t]){
      patch_above     <- which(population_K2[t, ] > K_temp_2[t]) # find which patch is above K
      patch_not_above <- which(!(population_K2[t, ] > K_temp_2[t])) # find patch not above K
      spillover       <- population_K2[t, patch_above] - K_temp_2[t] # set spillover to the difference between population and K
      
      population_K2[t, patch_above]     <- K_temp_2[t] # force patch above K to equal K
      population_K2[t, patch_not_above] <- population_K2[t, patch_not_above] + spillover # set the other patch equal to population plus spillover
      
      if(population_K2[t, patch_not_above] > K_temp_2[t]){ # need to check and make sure spillover to the other patch does not push the population over carrying capactiy
        population_K2[t, patch_not_above] <- K_temp_2[t] # if it does then set that patch to carrying capacity after spillover
      }
    }
    
  }
  
  outcome_population[iter, ] <- colMeans(population[(t-19):t,])
  outcome_harvest[iter, ] <- colMeans(harvest[(t-19):t,])
  
  outcome_population_r1[iter, ] <- colMeans(population_r1[(t-19):t,])
  outcome_harvest_r1[iter, ] <- colMeans(harvest_r1[(t-19):t,])
  
  outcome_population_r2[iter, ] <- colMeans(population_r2[(t-19):t,])
  outcome_harvest_r2[iter, ] <- colMeans(harvest_r2[(t-19):t,])
  
  outcome_population_r3[iter, ] <- colMeans(population_r3[(t-19):t,])
  outcome_harvest_r3[iter, ] <- colMeans(harvest_r3[(t-19):t,])
  
  outcome_population_K1[iter, ] <- colMeans(population_K1[(t-19):t,])
  outcome_harvest_K1[iter, ] <- colMeans(harvest_K1[(t-19):t,])
  
  outcome_population_K2[iter, ] <- colMeans(population_K2[(t-19):t,])
  outcome_harvest_K2[iter, ] <- colMeans(harvest_K2[(t-19):t,])
  
}



colnames(outcome_population) <- c("open_fish_notemp", "mpa_fish_notemp")
outcome_population <- as.data.frame(outcome_population)
colnames(outcome_harvest) <- c("open_harvest_notemp", "mpa_harvest_notemp")

colnames(outcome_population_r1) <- c("open_fish_r1", "mpa_fish_r1")
outcome_population_r1 <- as.data.frame(outcome_population_r1)
colnames(outcome_harvest_r1) <- c("open_harvest_r1", "mpa_harvest_r1")

colnames(outcome_population_r2) <- c("open_fish_r2", "mpa_fish_r2")
outcome_population_r2 <- as.data.frame(outcome_population_r2)
colnames(outcome_harvest_r2) <- c("open_harvest_r2", "mpa_harvest_r2")

colnames(outcome_population_r3) <- c("open_fish_r3", "mpa_fish_r3")
outcome_population_r3 <- as.data.frame(outcome_population_r3)
colnames(outcome_harvest_r3) <- c("open_harvest_r3", "mpa_harvest_r3")

colnames(outcome_population_K1) <- c("open_fish_K1", "mpa_fish_K1")
outcome_population_K1 <- as.data.frame(outcome_population_K1)
colnames(outcome_harvest_K1) <- c("open_harvest_K1", "mpa_harvest_K1")

colnames(outcome_population_K2) <- c("open_fish_K2", "mpa_fish_K2")
outcome_population_K2 <- as.data.frame(outcome_population_K2)
colnames(outcome_harvest_K2) <- c("open_harvest_K2", "mpa_harvest_K2")

#Fish population dataframe
outcome <- cbind(parameter_grid, outcome_population, outcome_population_r1, outcome_population_r2, outcome_population_r3,
                 outcome_population_K1, outcome_population_K2)

outcome$area_open <- map_dbl(outcome$patch_area, 1)
outcome$area_mpa <- map_dbl(outcome$patch_area, 2)

#Calculate weighted averages for each model version
outcome$wtavg_fish_notemp <- ((outcome$open_fish_notemp * outcome$area_open) + 
                                (outcome$mpa_fish_notemp * outcome$area_mpa)) / (outcome$area_open + outcome$area_mpa)

outcome$wtavg_fish_r1 <- ((outcome$open_fish_r1 * outcome$area_open) + 
                            (outcome$mpa_fish_r1 * outcome$area_mpa)) / (outcome$area_open + outcome$area_mpa)

outcome$wtavg_fish_r2 <- ((outcome$open_fish_r2 * outcome$area_open) + 
                            (outcome$mpa_fish_r2 * outcome$area_mpa)) / (outcome$area_open + outcome$area_mpa)

outcome$wtavg_fish_r3 <- ((outcome$open_fish_r3 * outcome$area_open) + 
                            (outcome$mpa_fish_r3 * outcome$area_mpa)) / (outcome$area_open + outcome$area_mpa)

outcome$wtavg_fish_K1 <- ((outcome$open_fish_K1 * outcome$area_open) + 
                            (outcome$mpa_fish_K1 * outcome$area_mpa)) / (outcome$area_open + outcome$area_mpa)

outcome$wtavg_fish_K2 <- ((outcome$open_fish_K2 * outcome$area_open) + 
                            (outcome$mpa_fish_K2 * outcome$area_mpa)) / (outcome$area_open + outcome$area_mpa)


outcome<-outcome[,c(1:14, 17:22, 15, 16)]

outcome_long <- pivot_longer(outcome, open_fish_notemp:wtavg_fish_K2, names_to = "model_version",
                             values_to = "population")

#total Fish - wt avg both patches combined
outcome_long_wtavg <- outcome_long %>%
  filter(model_version == "wtavg_fish_notemp" |
           model_version == "wtavg_fish_r1" |
           model_version == "wtavg_fish_r2" |
           model_version == "wtavg_fish_r3" |
           model_version == "wtavg_fish_K1" |
           model_version == "wtavg_fish_K2"
  )

outcome_long_wtavg$model_version <- factor(outcome_long_wtavg$model_version, 
                                           levels = c("wtavg_fish_notemp", "wtavg_fish_r1", "wtavg_fish_r2", 
                                                      "wtavg_fish_r3", "wtavg_fish_K1", "wtavg_fish_K2"),
                                           labels = c("No Temp", "r1", "r2", 
                                                      "r3", "K1", "K2"))

ggplot(outcome_long_wtavg, aes(x = area_mpa, y = population, col = model_version)) +
  geom_line() +
  facet_wrap(~S) +
  scale_color_viridis_d(name="Model version") +
  labs(x = "MPA area (proportion)", y = bquote("Wt avg fish density"~(g/m^2))) +
  theme_minimal() +
  ggtitle("Site fidelity") +
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1"))

#patch 1 fish only
outcome_long_open <- outcome_long %>%
  filter(model_version == "open_fish_notemp" |
           model_version == "open_fish_r1" |
           model_version == "open_fish_r2" |
           model_version == "open_fish_r3" |
           model_version == "open_fish_K1" |
           model_version == "open_fish_K2"
  )

outcome_long_open$model_version <- factor(outcome_long_open$model_version, 
                                          levels = c("open_fish_notemp", "open_fish_r1", "open_fish_r2", 
                                                     "open_fish_r3", "open_fish_K1", "open_fish_K2"),
                                          labels = c("No Temp", "r1", "r2", 
                                                     "r3", "K1", "K2"))

ggplot(outcome_long_open, aes(x = area_mpa, y = population, col = model_version)) +
  geom_line() +
  facet_wrap(~S, scales="free") +
  scale_color_viridis_d(name = "Model version") +
  labs(x = "MPA area (proportion)", y = bquote("Open area fish density"~(g/m^2))) +
  theme_minimal() +
  ggtitle("Site fidelity") +
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1"))


#patch 2 only
outcome_long_mpa <- outcome_long %>%
  filter(model_version == "mpa_fish_notemp" |
           model_version == "mpa_fish_r1" |
           model_version == "mpa_fish_r2" |
           model_version == "mpa_fish_r3" |
           model_version == "mpa_fish_K1" |
           model_version == "mpa_fish_K2"
  )

outcome_long_mpa$model_version <- factor(outcome_long_mpa$model_version, 
                                         levels = c("mpa_fish_notemp", "mpa_fish_r1", "mpa_fish_r2", 
                                                    "mpa_fish_r3", "mpa_fish_K1", "mpa_fish_K2"),
                                         labels = c("No Temp", "r1", "r2", 
                                                    "r3", "K1", "K2"))

ggplot(outcome_long_mpa, aes(x = area_mpa, y = population, col = model_version)) +
  geom_line() +
  facet_wrap(~S) +
  scale_color_viridis_d(name="Model version") +
  labs(x = "MPA area (proportion)", y = bquote("MPA fish density"~(g/m^2))) +
  theme_minimal() +
  #ggtitle("Site fidelity")
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1"))


ggarrange(p1+rremove("xlab"), p2, nrow=2, ncol=1, common.legend = TRUE, legend = "bottom")



#Harvest dataframe
outcome_harvest <- cbind(parameter_grid, outcome_harvest, outcome_harvest_r1, outcome_harvest_r2, outcome_harvest_r3,
                         outcome_harvest_K1, outcome_harvest_K2)

outcome_harvest$area_open <- map_dbl(outcome_harvest$patch_area, 1)
outcome_harvest$area_mpa <- map_dbl(outcome_harvest$patch_area, 2)
#outcome_harvest$fishing_p1 <- map_dbl(outcome_harvest$fishing_effort, 1)
#outcome_harvest$fishing_p2 <- map_dbl(outcome_harvest$fishing_effort, 2)

outcome_harvest_long <- pivot_longer(outcome_harvest, open_harvest_notemp:mpa_harvest_K2, names_to = "model_version",
                                     values_to = "harvest")


#patch 1 fish only
outcome_harvest_long_open <- outcome_harvest_long %>%
  filter(model_version == "open_harvest_notemp" |
           model_version == "open_harvest_r1" |
           model_version == "open_harvest_r2" |
           model_version == "open_harvest_r3" |
           model_version == "open_harvest_K1" |
           model_version == "open_harvest_K2"
  )

outcome_harvest_long_open$model_version <- factor(outcome_harvest_long_open$model_version, 
                                                  levels = c("open_harvest_notemp", "open_harvest_r1", "open_harvest_r2", 
                                                             "open_harvest_r3", "open_harvest_K1", "open_harvest_K2"),
                                                  labels = c("No Temp", "r1", "r2", 
                                                             "r3", "K1", "K2"))

ggplot(outcome_harvest_long_open, aes(x = area_mpa, y = harvest, col = model_version)) +
  geom_line() +
  facet_wrap(~S) +
  scale_color_viridis_d(name="Model version") +
  labs(x = "MPA area (proportion)", y = bquote("Harvest density"~(g/m^2))) +
  theme_minimal() +
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5), 
        legend.position = "bottom") +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1"))

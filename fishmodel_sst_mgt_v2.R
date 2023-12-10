library(tidyverse)
library(ggpubr)

source('dispersal.R')

projections <- read.csv("anomaly_df.csv")
SST_dev <- projections[['anomaly']]

timesteps <- length(SST_dev)
r = 0.3
K = c(101.3, 101.3)
S = 0.5

patch_area_sequences <- list(seq(0, 1, by = 0.1))
patch_area_grid <- do.call(expand.grid, patch_area_sequences)
patch_area_grid$Var2 <- 1 - patch_area_grid$Var1
patch_area_list <- split(patch_area_grid, 1:nrow(patch_area_grid))
number_patches <- ncol(patch_area_grid)

fishing_effort_sequences <- list(seq(0, 1, by = 0.1), 0)
fishing_effort_grid <- do.call(expand.grid, fishing_effort_sequences)
fishing_effort_list <- split(fishing_effort_grid, 1:nrow(fishing_effort_grid))
catchability <- 1

parameter_grid <- expand.grid(patch_area = patch_area_list,
                              fishing_effort = fishing_effort_list)

#saving outputs
outcome_population <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))
outcome_harvest_notemp <- matrix(0, ncol = 2, nrow = nrow(parameter_grid))

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
  tmp <- dispersal(as.numeric(parameter_grid[['patch_area']][[iter]]), number_patches, S)
  
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
    harvest[t,] <- population[t,] * (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population[t,] <- population[t,] * (1 - (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)))
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
    harvest_r1[t,] <- population_r1[t,] * (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population_r1[t,] <- population_r1[t,] * (1 - (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)))
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
    r_temp_2[t] <- 0.3 + 0*(SST_dev[t] - 1) + -0.0037*(SST_dev[t] - 1 )^2
    
    if(as.numeric(parameter_grid[['patch_area']][[iter]][1]) == 0 | as.numeric(parameter_grid[['patch_area']][[iter]][2]) == 0){
      patch_zero     <- which(as.numeric(parameter_grid[['patch_area']][[iter]]) == 0) # find which patch is above K
      population_r2[t-1, patch_zero]     <- 0 # force patch above K to equal K
    }
    population_r2[t,] <- population_r2[t-1,] + r_temp_2[t] * population_r2[t-1,] * (1 - population_r2[t-1,] / K)
    harvest_r2[t,] <- population_r2[t,] * (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population_r2[t,] <- population_r2[t,] * (1 - (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)))
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
    r_temp_3[t] <- 0.3 + 0*(SST_dev[t] + 1) + -0.0037*(SST_dev[t] + 1 )^2
    
    if(as.numeric(parameter_grid[['patch_area']][[iter]][1]) == 0 | as.numeric(parameter_grid[['patch_area']][[iter]][2]) == 0){
      patch_zero     <- which(as.numeric(parameter_grid[['patch_area']][[iter]]) == 0) 
      population_r3[t-1, patch_zero]     <- 0
    }
    population_r3[t,] <- population_r3[t-1,] + r_temp_3[t] * population_r3[t-1,] * (1 - population_r3[t-1,] / K)
    harvest_r3[t,] <- population_r3[t,] * (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population_r3[t,] <- population_r3[t,] * (1 - (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)))
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
    harvest_K1[t,] <- population_K1[t,] * (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population_K1[t,] <- population_K1[t,] * (1 - (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)))
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
    if(as.numeric(parameter_grid[['patch_area']][[iter]][,1]) == 0){
      population_K1[t,1] <- 0
    } else if (as.numeric(parameter_grid[['patch_area']][[iter]][,2]) == 0) {
      population_K1[t,2] <- 0
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
    harvest_K2[t,] <- population_K2[t,] * (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)) * 
      as.numeric(parameter_grid[['patch_area']][[iter]])
    population_K2[t,] <- population_K2[t,] * (1 - (1 - exp(-as.numeric(parameter_grid[['fishing_effort']][[iter]]) * catchability)))
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
    if(as.numeric(parameter_grid[['patch_area']][[iter]][,1]) == 0){
      population_K2[t,1] <- 0
    } else if (as.numeric(parameter_grid[['patch_area']][[iter]][,2]) == 0) {
      population_K2[t,2] <- 0
    }
  }
  
  outcome_population[iter, ] <- colMeans(population[(t-19):t,])
  outcome_harvest_notemp[iter, ] <- colMeans(harvest[(t-19):t,])
  
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
colnames(outcome_harvest_notemp) <- c("open_harvest_notemp", "mpa_harvest_notemp")

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
outcome$fishing_p1 <- map_dbl(outcome$fishing_effort, 1)
outcome$fishing_p2 <- map_dbl(outcome$fishing_effort, 2)


###########.   MPA FISH BIOMASS.  #####################

outcome_mpa <- outcome %>%
  select(mpa_fish_notemp, mpa_fish_r1, mpa_fish_r2, mpa_fish_r3, mpa_fish_K1, mpa_fish_K2, area_open, area_mpa, fishing_p1)

outcome_mpa_long <- outcome_mpa %>%
  pivot_longer(cols = mpa_fish_notemp: mpa_fish_K2, names_to = "model_version", values_to = "mpafish_biomass")

outcome_mpa_long$model_version <- factor(outcome_mpa_long$model_version, 
                                         levels = c("mpa_fish_notemp", "mpa_fish_r1", "mpa_fish_r2", 
                                                    "mpa_fish_r3", "mpa_fish_K1", "mpa_fish_K2"),
                                         labels = c("Baseline", "r1", "r2", 
                                                    "r3", "K1", "K2"))
p1<-ggplot(outcome_mpa_long %>%
         filter(fishing_p1 == 0.1 | fishing_p1 == 0.5 | fishing_p1 == 0.9), 
       aes(x = area_mpa, y = mpafish_biomass, col = model_version)) +
  geom_line(lwd=0.75) +
  facet_wrap(~fishing_p1) +
  scale_color_viridis_d(name="Model version") +
  labs(x = "MPA area (proportion)", y = bquote("MPA fish biomass"~(g/m^2))) +
  theme_minimal() +
  ggtitle("A.") +
  theme(text = element_text(size=20), legend.position = "bottom") +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1"))



#RELATIVE BIOMASS

mpafish_sum <- outcome_mpa %>% filter(area_open == 0) %>%
  group_by(fishing_p1) %>%
  summarize(mpa_fish_notemp, mpa_fish_r1, mpa_fish_r2, mpa_fish_r3, mpa_fish_K1, mpa_fish_K2)

condition_mapping <- data.frame(
  fishing_p1 = unique(outcome_mpa$fishing_p1),
  mpafish_full_mpa = c(mpafish_sum[1:11, 2]), 
  mpafish_full_mpa_r1 = c(mpafish_sum[1:11, 3]),
  mpafish_full_mpa_r2 = c(mpafish_sum[1:11, 4]), 
  mpafish_full_mpa_r3 = c(mpafish_sum[1:11, 5]), 
  mpafish_full_mpa_k1 = c(mpafish_sum[1:11, 6]), 
  mpafish_full_mpa_k2 = c(mpafish_sum[1:11, 7])
)

colnames(condition_mapping) <- c("fishing_p1", "mpafish_full_mpa" , "mpafish_full_mpa_r1", "mpafish_full_mpa_r2", 
                                 "mpafish_full_mpa_r3", "mpafish_full_mpa_K1", "mpafish_full_mpa_K2")

outcome_mpa_final <- merge(outcome_mpa, condition_mapping, by = "fishing_p1")

outcome_mpa_final$rel_mpafish_no_temp <- replace(outcome_mpa_final$mpa_fish_notemp / outcome_mpa_final$mpafish_full_mpa, is.nan(outcome_mpa_final$mpa_fish_notemp / outcome_mpa_final$mpafish_full_mpa), 0)
outcome_mpa_final$rel_mpafish_r1 <- replace(outcome_mpa_final$mpa_fish_r1 / outcome_mpa_final$mpafish_full_mpa_r1, is.nan(outcome_mpa_final$mpa_fish_r1 / outcome_mpa_final$mpafish_full_mpa_r1), 0)
outcome_mpa_final$rel_mpafish_r2 <- replace(outcome_mpa_final$mpa_fish_r2 / outcome_mpa_final$mpafish_full_mpa_r2, is.nan(outcome_mpa_final$mpa_fish_r2 / outcome_mpa_final$mpafish_full_mpa_r2), 0)
outcome_mpa_final$rel_mpafish_r3 <- replace(outcome_mpa_final$mpa_fish_r3 / outcome_mpa_final$mpafish_full_mpa_r3, is.nan(outcome_mpa_final$mpa_fish_r3 / outcome_mpa_final$mpafish_full_mpa_r3), 0)
outcome_mpa_final$rel_mpafish_K1 <- replace(outcome_mpa_final$mpa_fish_K1 / outcome_mpa_final$mpafish_full_mpa_K1, is.nan(outcome_mpa_final$mpa_fish_K1 / outcome_mpa_final$mpafish_full_mpa_K1), 0)
outcome_mpa_final$rel_mpafish_K2 <- replace(outcome_mpa_final$mpa_fish_K2 / outcome_mpa_final$mpafish_full_mpa_K2, is.nan(outcome_mpa_final$mpa_fish_K2 / outcome_mpa_final$mpafish_full_mpa_K2), 0)


outcome_mpa_final_long <- outcome_mpa_final %>%
  select(fishing_p1, area_open, area_mpa, rel_mpafish_no_temp:rel_mpafish_K2) %>%
  pivot_longer(cols = rel_mpafish_no_temp:rel_mpafish_K2, names_to = "model_version",
               values_to = "rel_mpa_biomass")

outcome_mpa_final_long$model_version <- factor(outcome_mpa_final_long$model_version, 
                                           levels = c("rel_mpafish_no_temp", "rel_mpafish_r1", "rel_mpafish_r2", 
                                                      "rel_mpafish_r3", "rel_mpafish_K1", "rel_mpafish_K2"),
                                           labels = c("Baseline", "r1", "r2", 
                                                      "r3", "K1", "K2"))

p2<-ggplot(outcome_mpa_final_long %>% filter(fishing_p1 == 0.1 |
                                       fishing_p1 == 0.5 |
                                       fishing_p1 == 0.9), 
       aes(x=area_mpa, y=rel_mpa_biomass, col=model_version)) +
  facet_wrap(~fishing_p1, scales="free") +
  geom_line(lwd=1) +
  scale_color_viridis_d(name="Model version") +
  labs(x = "MPA area (proportion)", y = "Relative MPA fish biomass",
       caption = "relative to MPA =1") +
  theme_minimal() +
  ggtitle("B.") +
  theme(text = element_text(size=20), legend.position = "bottom", plot.caption = element_text(hjust = 0)) +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1"))


figure<-ggarrange(p1+rremove("xlab"), p2+rremove("xlab"), 
                  nrow=2, ncol=1, common.legend = TRUE, legend = "top")


annotate_figure(figure, bottom = text_grob("MPA area (proportion)",
                                           size = 20))





#zoomed in
ggplot(outcome_mpa_final_long %>% filter(fishing_p1 == 0.1 |
                                           fishing_p1 == 0.5 |
                                           fishing_p1 == 0.9) %>%
         filter(area_mpa > 0.5),
       aes(x=area_mpa, y=rel_mpa_biomass, col=model_version)) +
  facet_wrap(~fishing_p1, scales="free") +
  geom_line() +
  scale_color_viridis_d(name="Model version") +
  labs(x = "MPA area (proportion)", y = "MPA fish biomass relative to MPA = 1") +
  theme_minimal() +
  ggtitle("Fishing effort") +
  theme(text = element_text(size=20), plot.title = element_text(hjust = 0.5), legend.position = "bottom") 

#+
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1"))


  
  
  
  
  

##############HARVEST############################
  
outcome_harvest <- cbind(parameter_grid, outcome_harvest_notemp, outcome_harvest_r1, outcome_harvest_r2, outcome_harvest_r3,
                 outcome_harvest_K1, outcome_harvest_K2)

outcome_harvest$area_open <- map_dbl(outcome_harvest$patch_area, 1)
outcome_harvest$area_mpa <- map_dbl(outcome_harvest$patch_area, 2)
outcome_harvest$fishing_p1 <- map_dbl(outcome_harvest$fishing_effort, 1)
outcome_harvest$fishing_p2 <- map_dbl(outcome_harvest$fishing_effort, 2)


outcome_harv2 <- outcome_harvest %>%
  select(open_harvest_notemp, open_harvest_r1, open_harvest_r2, open_harvest_r3, 
         open_harvest_K1, open_harvest_K2,
         area_open, area_mpa, fishing_p1)


outcome_harv_long <- outcome_harv2 %>%
  pivot_longer(cols = open_harvest_notemp: open_harvest_K2, names_to = "model_version", values_to = "harvest")

outcome_harv_long$model_version <- factor(outcome_harv_long$model_version, 
                                         levels = c("open_harvest_notemp", "open_harvest_r1", "open_harvest_r2", 
                                                    "open_harvest_r3", "open_harvest_K1", "open_harvest_K2"),
                                         labels = c("Baseline", "r1", "r2", 
                                                    "r3", "K1", "K2"))

p1<-ggplot(outcome_harv_long %>%
         filter(fishing_p1 == 0.1 | fishing_p1 == 0.5 | fishing_p1 == 0.9), 
       aes(x = area_mpa, y = harvest, col = model_version)) +
  geom_line(lwd=0.75) +
  facet_wrap(~fishing_p1) +
  scale_color_viridis_d(name="Model version") +
  labs(x = "MPA area (proportion)", y = bquote("Harvest"~(g/m^2))) +
  theme_minimal() +
  ggtitle("A.") +
  theme(text = element_text(size=20),
        legend.position = "bottom") +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1"))



##RELATIVE HARVEST
#Indexing to create havest relative to harvest at MPA = 0

harv_sum <- outcome_harvest %>% filter(area_open == 1.0) %>%
  group_by(fishing_p1) %>%
  summarize(open_harvest_notemp, open_harvest_r1, open_harvest_r2, open_harvest_r3, open_harvest_K1, open_harvest_K2)
  
condition_mapping <- data.frame(
  fishing_p1 = unique(outcome_harvest$fishing_p1),
  harv_no_mpa = c(harv_sum[1:11, 2]), 
  harv_no_mpa_r1 = c(harv_sum[1:11, 3]),
  harv_no_mpa_r2 = c(harv_sum[1:11, 4]), 
  harv_no_mpa_r3 = c(harv_sum[1:11, 5]), 
  harv_no_mpa_k1 = c(harv_sum[1:11, 6]), 
  harv_no_mpa_k2 = c(harv_sum[1:11, 7])
)

colnames(condition_mapping) <- c("fishing_p1", "harv_no_mpa" , "harv_no_mpa_r1", "harv_no_mpa_r2", 
                                 "harv_no_mpa_r3", "harv_no_mpa_K1", "harv_no_mpa_K2")

outcome_final <- merge(outcome_harvest, condition_mapping, by = "fishing_p1")

outcome_final$rel_harv_no_temp <- outcome_final$open_harvest_notemp / outcome_final$harv_no_mpa
outcome_final$rel_harv_r1 <- outcome_final$open_harvest_r1 / outcome_final$harv_no_mpa_r1
outcome_final$rel_harv_r2 <- outcome_final$open_harvest_r2 / outcome_final$harv_no_mpa_r2
outcome_final$rel_harv_r3 <- outcome_final$open_harvest_r3 / outcome_final$harv_no_mpa_r3
outcome_final$rel_harv_K1 <- outcome_final$open_harvest_K1 / outcome_final$harv_no_mpa_K1
outcome_final$rel_harv_K2 <- outcome_final$open_harvest_K2 / outcome_final$harv_no_mpa_K2


outcome_final$rel_harv_no_temp[is.nan(outcome_final$rel_harv_no_temp)] <- 0
outcome_final$rel_harv_r1[is.nan(outcome_final$rel_harv_r1)] <- 0
outcome_final$rel_harv_r2[is.nan(outcome_final$rel_harv_r2)] <- 0
outcome_final$rel_harv_r3[is.nan(outcome_final$rel_harv_r3)] <- 0
outcome_final$rel_harv_K1[is.nan(outcome_final$rel_harv_K1)] <- 0
outcome_final$rel_harv_K1[is.nan(outcome_final$rel_harv_K1)] <- 0


outcome_final_long <- outcome_final %>%
  select(fishing_p1, area_open, area_mpa, rel_harv_no_temp:rel_harv_K2) %>%
  pivot_longer(cols = rel_harv_no_temp:rel_harv_K2, names_to = "model_version",
                                   values_to = "rel_harvest")

outcome_final_long$model_version <- factor(outcome_final_long$model_version, 
                                           levels = c("rel_harv_no_temp", "rel_harv_r1", "rel_harv_r2", 
                                                      "rel_harv_r3", "rel_harv_K1", "rel_harv_K2"),
                                           labels = c("Baseline", "r1", "r2", 
                                                      "r3", "K1", "K2"))

p2<-ggplot(outcome_final_long %>% filter(fishing_p1 == 0.1 |
                                       fishing_p1 == 0.5 |
                                       fishing_p1 == 0.9
                                       ), aes(x=area_mpa, y=rel_harvest, col=model_version)) +
  facet_wrap(~fishing_p1, scales="free") +
  geom_line(lwd=0.75) +
  scale_color_viridis_d(name="Model version") +
  labs(x = "MPA area (proportion)", y = "Harvest relative to MPA = 0") +
  theme_minimal() +
  ggtitle("B.") +
  theme(text = element_text(size=20), legend.position = "bottom") +
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1"))



figure<-ggarrange(p1+rremove("xlab"), p2+rremove("xlab"), 
                  nrow=2, ncol=1, common.legend = TRUE, legend = "top")

annotate_figure(figure, bottom = text_grob("MPA area (proportion)",
                                           size = 20))


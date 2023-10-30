#Schaefer model
patch_area <- c(1, 0)
number_patches <- length(patch_area)
timesteps <- 80
r = 0.3
K = 101.3

projections <- read.csv("projections.csv")
SST_dev <- projections[['anomaly']]
SST_dev_higheropttemp <- projections[['anomaly_higheropttemp']]
SST_dev_loweropttemp <- projections[['anomaly_loweropttemp']]

mean_temps_crop <- read.csv("mean_temps_crop.csv")
mean_temps_series <- mean_temps_crop[['extracted_mean_temps']]

population <- array(NA, dim = c(timesteps, number_patches))
population[1,] <- c(10,0)

population_r1 <- array(NA, dim = c(timesteps, number_patches))
population_r1[1,] <- c(10,0)

population_r2 <- array(NA, dim = c(timesteps, number_patches))
population_r2[1,] <- c(10,0)

population_r3 <- array(NA, dim = c(timesteps, number_patches))
population_r3[1,] <- c(10,0)

population_K1 <- array(NA, dim = c(timesteps, number_patches))
population_K1[1,] <- c(10,0)

population_K2 <- array(NA, dim = c(timesteps, number_patches))
population_K2[1,] <- c(10,0)


r_temp_1 <- array(NA, dim = c(timesteps, 1))
r_temp_2 <- array(NA, dim = c(timesteps, 1))
r_temp_3 <- array(NA, dim = c(timesteps, 1))

K_temp_1 <- array(NA, dim = c(timesteps, 1))
K_temp_2 <- array(NA, dim = c(timesteps, 1))




for (t in 2:timesteps) {
  for (i in 1:number_patches){
    
    #no temp
    population[t,i] <- population[t-1, i] + r * population[t-1, i] * (1 - population[t-1,i] / K)
    
    #temp dependent r - "current" temp as optimal
    r_temp_1[t] <- 0.3 + 0*SST_dev[t] + -0.0037*SST_dev[t]^2
    population_r1[t,i] <- population_r1[t-1, i] + r_temp_1[t] * population_r1[t-1, i] * (1 - population_r1[t-1,i] / K)
    
    #temp dependent r - higher optimal temp
    r_temp_2[t] <- a + b*SST_dev_higheropttemp[t] + c*SST_dev_higheropttemp[t]^2
    population_r2[t,i] <- population_r2[t-1, i] + r_temp_2[t] * population_r2[t-1, i] * (1 - population_r2[t-1,i] / K)
    
    #temp dependent r - lower optimal temp
    r_temp_3[t] <- a + b*SST_dev_loweropttemp[t] + c*SST_dev_loweropttemp[t]^2
    population_r3[t,i] <- population_r3[t-1, i] + r_temp_3[t] * population_r3[t-1, i] * (1 - population_r3[t-1,i] / K)
    
    #temp dependent K - piece wise linear function
    K_temp_1[t] <- -2.1408 * mean_temps_series[t] + 160.253358
    if (K_temp_1[t] < 10) {
      K_temp_1[t] <- 10
    } else if (K_temp_1[t] > 101.3) {
      K_temp_1[t] <- 101.3
    }
    population_K1[t,i] <- population_K1[t-1, i] + r * population_K1[t-1, i] * (1 - population_K1[t-1,i] / K_temp_1[t])
    
    #temp dependent K - quadratic function
    K_temp_2[t] <- 101.3 + 0*SST_dev[t] + -0.7*SST_dev[t]^2
    if (K_temp_2[t] < 10) {
      K_temp_2[t] <- 10 }
    else if (K_temp_2[t] > 101.3) {
      K_temp_2[t] <- 101.3
    }
    population_K2[t,i] <- population_K2[t-1, i] + r * population_K2[t-1, i] * (1 - population_K2[t-1,i] / K_temp_2[t])
  }
}



plot(population[,1])
plot(population_r1[,1])
plot(population_r2[,1])
plot(population_r3[,1])
plot(population_K1[,1])
plot(population_K2[,1])

time <- c(1:80)
df <-cbind(NoTemp = population[,1], r1 = population_r1[,1], r2 = population_r2[,1], r3 = population_r3[,1],
           K1 = population_K1[,1], K2 = population_K2[,1], time)

df_long <- pivot_longer(as.data.frame(df), NoTemp:K2, names_to = "Temp_version", values_to = "Population")

df_long$Temp_version <- factor(df_long$Temp_version, levels = c("NoTemp", "r1", "r2", "r3", "K1", "K2"),
                               labels = c("No Temp", "r1 (optimal temp: current)", "r2 (higher optimal temp)", 
                                          "r3 (lower optimal temp)", "K1 (linear)", "K2 (quadratic)"))

ggplot(df_long, aes(x=time, y=Population, col=Temp_version)) +
  geom_line() +
  scale_color_viridis_d()+
  theme_minimal()

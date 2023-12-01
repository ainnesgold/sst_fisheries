library(tidyverse)
library(ggpubr)
library(grid)

projections <- read.csv("anomaly_df.csv")
SST_dev <- projections[['anomaly']]
#SST_dev_loweropttemp <- projections[['anomaly']] + 1 #larger anomaly = lower optimal temp
#SST_dev_higheropttemp <- projections[['anomaly']] - 1 #smaller anomaly = higher optimal temp

patch_area <- c(1, 0)
number_patches <- length(patch_area)
timesteps <- length(SST_dev)
r = 0.3
K = 101.3

#mean_temps_crop <- read.csv("mean_temps_crop.csv")
#mean_temps_series <- mean_temps_crop[['extracted_mean_temps']]

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
    
    #no temp
    population[t,] <- population[t-1,] + r * population[t-1,] * (1 - population[t-1,] / K)
    
    #temp dependent r - "current" temp as optimal
    r_temp_1[t] <- 0.3 + 0*SST_dev[t] + -0.0037*SST_dev[t]^2
    population_r1[t,] <- population_r1[t-1,] + r_temp_1[t] * population_r1[t-1,] * (1 - population_r1[t-1,] / K)
    
    #temp dependent r - higher optimal temp
    r_temp_2[t] <- 0.3 + 0*(SST_dev[t] - 1) + -0.0037*(SST_dev[t] - 1 )^2
    population_r2[t,] <- population_r2[t-1,] + r_temp_2[t] * population_r2[t-1,] * (1 - population_r2[t-1,] / K)
    
    #temp dependent r - lower optimal temp
    r_temp_3[t] <- 0.3 + 0*(SST_dev[t] + 1) + -0.0037*(SST_dev[t] + 1 )^2
    population_r3[t,] <- population_r3[t-1,] + r_temp_3[t] * population_r3[t-1,] * (1 - population_r3[t-1,] / K)
    
    #temp dependent K - piece wise linear function
    K_temp_1[t] <- -4.95243768 * SST_dev[t] + 101.3
    if (K_temp_1[t] < 10) {
      K_temp_1[t] <- 10
    } else if (K_temp_1[t] > 101.3) {
      K_temp_1[t] <- 101.3
    }
    population_K1[t,] <- population_K1[t-1,] + r * population_K1[t-1,] * (1 - population_K1[t-1,] / K_temp_1[t])
    
    #temp dependent K - quadratic function
    K_temp_2[t] <- 101.3 + 0*SST_dev[t] + -0.7*SST_dev[t]^2
    if (K_temp_2[t] < 10) {
      K_temp_2[t] <- 10 }
    else if (K_temp_2[t] > 101.3) {
      K_temp_2[t] <- 101.3
    }
    population_K2[t,] <- population_K2[t-1,] + r * population_K2[t-1,] * (1 - population_K2[t-1,] / K_temp_2[t])
}



plot(population[,1])
plot(population_r1[,1])
plot(population_r2[,1])
plot(population_r3[,1])
plot(population_K1[,1])
plot(population_K2[,1])

time <- c(1:length(SST_dev))
df <-cbind(NoTemp = population[,1], r1 = population_r1[,1], r2 = population_r2[,1], r3 = population_r3[,1],
           K1 = population_K1[,1], K2 = population_K2[,1], time)

df_long <- pivot_longer(as.data.frame(df), NoTemp:K2, names_to = "Temp_version", values_to = "Population")

df_long$Temp_version <- factor(df_long$Temp_version, levels = c("NoTemp", "r1", "r2", "r3", "K1", "K2"),
                               labels = c("Baseline", "r1", "r2", 
                                          "r3", "K1", "K2"))

p1<-ggplot(df_long, aes(x=time, y=Population, col=Temp_version)) +
  geom_line() +
  scale_color_viridis_d(name = "Model version")+
  theme_minimal() +
  ggtitle("A.") +
  theme(text = element_text(size=20), legend.position = "top") +
  labs(x = "Years", y = bquote("Fish biomass"~(g/m^2)))


#zoomed in
p2<-ggplot(df_long %>% filter(time > 15 & time < 25), aes(x=time, y=Population, col=Temp_version)) +
  geom_line(lwd=1, position=position_dodge(width=0.2)) +
  scale_color_viridis_d(name = "Model version")+
  theme_minimal() +
  ggtitle("B.") +
  theme(text = element_text(size=20), legend.position = "top") +
  labs(x = "Years", y = bquote("Fish biomass"~(g/m^2)))


figure<-ggarrange(p1  + rremove("xlab"),p2+rremove("ylab") + rremove("xlab"),nrow=1, common.legend = TRUE, legend = "top")

annotate_figure(figure, bottom = text_grob("Years",
                                  size = 20))



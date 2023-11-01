# Code from: https://rpubs.com/markpayne/358146
#Dataset was too large to upload to github but can be downloaded here
#CMIP6: monthly, single levels, SST, CESM2 (USA), whole temporal range, whole available region
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/projections-cmip6?tab=form

library(tidyverse)
library(raster)
library(ncdf4)

#2015-2064
##SSP2-4.5
fname <- "CMIP6 trials/4.5/tos_Omon_CESM2_ssp245_r4i1p1f1_gn_201501-206412_v20200528.nc"

temp_file <- brick(fname)
plot(temp_file[[1]])

lon.pts <- seq(202,202.99,by=0.01) #~long of Oahu edges (calculated as 360 - lon to get degrees east). 202-202.99 by 0.01
lat.pts <- seq(201.5, 202.49, by=0.01)  #degrees north. equator is 180 + oahu latittude edges (21) = 201, 201.99
extract.pts <- cbind(lon.pts,lat.pts)
ann.extract <- extract(temp_file,extract.pts,method="bilinear")
extracted_number_years <- (length(ann.extract) / 1200.0)

plot(crop(temp_file[[1]], extract.pts))


#2015 to 2064
years <- array(2015:2064)
extracted_mean_temps <- c()

for (t in 0:(floor(extracted_number_years)-1)) {
  year_step <- t * 100
  b <- ann.extract[(year_step+1):(year_step+100)]
  extracted_mean_temp <- mean(b)
  extracted_mean_temps <- append(extracted_mean_temps, extracted_mean_temp)
}

extracted_mean_temps<- cbind(years, extracted_mean_temps)
mean_temps_2064<-as.data.frame(extracted_mean_temps)

#2065-2100
fname2 <- "CMIP6 trials/4.5/tos_Omon_CESM2_ssp245_r4i1p1f1_gn_206501-210012_v20200528.nc"
temp_file2 <- brick(fname2)
ann.extract <- extract(temp_file2,extract.pts,method="bilinear")
extracted_number_years <- (length(ann.extract) / 1200.0)

years <- array(2065:2100)
extracted_mean_temps <- c()

for (t in 0:(floor(extracted_number_years)-1)) {
  year_step <- t * 100
  b <- ann.extract[(year_step+1):(year_step+100)]
  extracted_mean_temp <- mean(b)
  extracted_mean_temps <- append(extracted_mean_temps, extracted_mean_temp)
}

extracted_mean_temps<- cbind(years, extracted_mean_temps)
mean_temps_2100<-as.data.frame(extracted_mean_temps)


mean_temps<-rbind(mean_temps_2064, mean_temps_2100)
mean_temps <- as.data.frame(mean_temps)



#Calculate future years deviation from present
#historical data, too big for github but can be downloaded from:
#https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html HadISST netCDF at bottom of page
lon.pts <- seq(157.5, 158.5, by = 1)
lat.pts <- seq(21, 22, by = 1)
extract.pts <- cbind(lon.pts,lat.pts)

histname <- "CMIP6 trials/HadISST_sst.nc"
hist_file <- brick(histname)

##LOOP TO CALCULATE AVG TEMP FOR EACH YEAR 
ann.extract <- extract(hist_file,extract.pts,method="bilinear")
hist_number_years <- (length(ann.extract) / 24.0)

years <- array(1870:2020)
hist_mean_temps <- c()
hist_max_temps <- c()


for (t in 0:(floor(hist_number_years)-1)) {
  year_step <- t * 24
  b <- ann.extract[(year_step+1):(year_step+24)]
  hist_mean_temp <- mean(b)
  hist_mean_temps <- append(hist_mean_temps, hist_mean_temp)
  hist_max_temp <- max(b)
  hist_max_temps <- append(hist_max_temps, hist_max_temp)
}
hist_mean_temps<- cbind(years, hist_mean_temps)
hist_mean_temps<-as.data.frame(hist_mean_temps)

hist_max_temps<- cbind(years, hist_max_temps)
hist_max_temps<-as.data.frame(hist_max_temps)


current_temp <- hist_mean_temps %>%
  filter(years >= 1991,
         years <= 2021) %>%
  summarize(avgtemp = mean(hist_mean_temps))

current_max_temp <- hist_max_temps %>%
  filter(years >= 1991,
         years <= 2021) %>%
  summarize(avgmaxtemp = mean(hist_max_temps))

#create dataframe of future values only, starting at 2021. I want 50 years right now to match the model runs.
mean_temps_crop <- mean_temps %>%
  filter(years > 2020)

projections <- mean_temps %>%
  filter(years > 2020) %>%
  summarise(anomaly = extracted_mean_temps - current_temp$avgtemp,
            anomaly_higheropttemp = extracted_mean_temps - (current_temp$avgtemp + 1),
            anomaly_loweropttemp = extracted_mean_temps - (current_temp$avgtemp - 1))

write.csv(projections, "projections.csv")
write.csv(mean_temps_crop, "mean_temps_crop.csv")


rcp4.5<-mean_temps

#supplemental temp figure
ggplot(data=rcp4.5 %>% filter (years > 2019, years < 2071), aes(x=years, y=extracted_mean_temps)) +
  geom_point() +
  geom_line() +
  labs(x="Year", y="°C") +
  theme_minimal() +
  theme(text = element_text(size=20))

#get a 30 yr avg for 2040 - future
 rcp4.5 %>% filter(years >=2025 | years <= 2055) %>%
   summarise(mean(extracted_mean_temps))




##plotting historical mean and max temps
hist_mean_temps
hist_temps <- left_join(hist_mean_temps, hist_max_temps)


hist_long <- hist_temps %>% 
  pivot_longer(cols = "hist_mean_temps":"hist_max_temps", names_to = "measure", values_to = "value")

ggplot(data = hist_long %>%
         filter(years >= 1991, years <2021), 
       aes(x=years, y=value, color=measure)) +
  geom_point() + 
  geom_line() +
  scale_color_viridis_d(name = "", labels = c("max temperature", "mean temperature")) +
  theme_minimal() +
  labs(x="Year", y = "°C") +
  theme(legend.position = "top") +
  theme(text = element_text(size=20))



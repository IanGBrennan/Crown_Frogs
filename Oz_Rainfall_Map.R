library(rnaturalearth)
library(raster)
library(nasapower)
library(viridis)

# Following the method here: https://rspatialdata.github.io/rainfall.html

# Fetch Climate Data
# Using lonlat equivalent to Australia, we will set temporal_api = "CLIMATOLOGY". 
# Since the maximum processed area is 4.5 x 4.5 degrees (100 points), we will do this in step-by-step.

flag <- 1
for (i in seq(110, 150, 5)) {
  for (j in seq(-40, -10, 5)) {
    climate_avg_temp <- get_power(community = "AG",
                                  pars = "PRECTOTCORR",
                                  lonlat = c(i, (j - 5), (i + 5), j),
                                  temporal_api = "CLIMATOLOGY")
    if (flag == 1) {
      climate_avg <- climate_avg_temp
      flag <- 0
    } else{
      climate_avg <- rbind(climate_avg, climate_avg_temp)
    }
  }
}


climate_avg %>% datatable(extensions = c("Scroller", "FixedColumns"), options = list(
  deferRender = TRUE,
  scrollY = 350,
  scrollX = 350,
  dom = "t",
  scroller = TRUE,
  fixedColumns = list(leftColumns = 3)
))

# this gives you the mm/day for each month of the year
# we need to multiple each month by the number of days, then sum it to get annual rainfall

climate_avg$JAN <- climate_a

days.per.month <- c(31,28,31,30,31,30, 31, 31, 30, 31, 30, 31)

for (j in 1:length(days.per.month)){
  climate_avg[,3+j] <- climate_avg[,3+j]*days.per.month[j]
}
climate_avg$AVG_ANN <- rowSums(climate_avg[ , c(4:15)])
climate_avg$logAVG_ANN <- log(climate_avg$AVG_ANN)

# Getting world map
map <- ne_countries(country = 'australia', returnclass = "sf")

# Converting data to raster
r <- rasterFromXYZ(climate_avg[, c("LON", "LAT", "logAVG_ANN", "AVG_ANN")])

# Converting the raster into a data.frame
r_df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)

# Plot
ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = AVG_ANN)) +
  scale_fill_distiller(palette="Spectral", direction=1) +
  #geom_hex(data = r_df, aes(x=x, y=y, fill=logAVG_ANN), stat="identity") + 
  geom_sf(data = map, inherit.aes = FALSE, fill = NA) +
  #scale_fill_viridis(option="E", direction=-1) +
  labs(
    title = "Rainfall in inches",
    fill = "Annual Rainfall",
    subtitle = "Annual rainfall in Australia"
  ) + 
  labs(x = "Longitude", y = "Latitude")
  

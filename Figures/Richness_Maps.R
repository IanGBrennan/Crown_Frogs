# packages
library(galah)      # To download species data
#library(rayshader)  # For 3d rendering
library(tidyverse)  # Data wrangling
library(here)       # Safe paths
library(sf)         # Spatial features
library(ozmaps)     # For map of oz


# get a map and project to WGS84
oz_wgs84 <- ozmap_data(data = "country") |>
  st_transform(crs = st_crs("WGS84"))

## check map
ggplot(oz_wgs84) + geom_sf()


# create grid
oz_grid <- st_make_grid(oz_wgs84,
                        what = "polygons",
                        cellsize = 1.0,
                        square = FALSE,
                        flat_topped = TRUE)

# subset to grid cells that are within land
keep_hexes <- st_intersects(oz_grid, oz_wgs84)
keep_hexes <- as.data.frame(keep_hexes)$row.id
oz_grid <- oz_grid[keep_hexes]

## check
ggplot() +
  geom_sf(data = oz_wgs84) +
  geom_sf(data = oz_grid, fill = NA, color = "red")


get_counts <- function(hexagon){
  
  # convert to wkt
  wkt_string <- st_as_text(oz_grid[[hexagon]]) %>%
    sub(")))", "))", .) %>%
    sub("POLYGON ", "POLYGON", .)
  
  # get counts
  result <- galah_call() |>
    galah_geolocate(wkt_string) |>
    galah_identify("Myobatrachidae") |>
    galah_filter(profile = "ALA",
                 decimalLongitude > 110 & decimalLongitude < 160,
                 country == "Australia") |>
    atlas_counts(type = "species", # get species counts
                 limit = NULL)
  
  # light formatting to catch errors
  if(is.null(result)){
    tibble(count = NA, id = hexagon)
  }else{
    result$id <- hexagon
    result
  }
}

# download number of species for each polygon
counts_list <- map(seq_along(oz_grid), get_counts)

# bind lists to data frame
counts_df <- map_dfr(counts_list, rbind)


# convert to tibble, attach counts
oz_df <- st_as_sf(oz_grid)
oz_df$count <- NA
oz_df$count[counts_df$id] <- counts_df$count


# See top hexagons
oz_df %>%
  arrange(desc(count)) %>%
  head(10L)

#micro_oz_df <- oz_df
#hylid_oz_df <- oz_df
#myo_oz_df <- oz_df

# plot it all together
myo_map <- ggplot() +
  geom_sf(
    data = oz_df,
    mapping = aes(fill = count), # log10 + 1 transformed
    alpha = 1,
    color = NA) +
  scale_fill_distiller(name = "Number of Myobatrachid species",
                       #type = "seq",
                       direction = 1,
                       #limits = c(0,3),
                       #labels = c("10", "100"),
                       palette = "Blues")
                       # edit legend to be horizontal-bottom
                       #guide = guide_colorsteps(direction = "horizontal",
                       #                         label.position = "top",
                       #                         title.position = "bottom",
                       #                         title.hjust = 0.5)
  ) +
  # add map
  geom_sf(data = oz_wgs84,
          color = NA,
          fill = NA)  +
  # crop map
  coord_sf(xlim = c(110, 155), 
           ylim = c(-45, -10)) #+
  # Adjust text and make aesthetic more minimal
  #theme(title = element_text(face = "bold"),
  #      legend.title = element_text(size = 19),
  #      legend.position = "bottom",
  #      legend.key.width = unit(28, 'mm'),
  #      legend.text = element_text(size = 16),
  #      plot.background = element_rect(fill = 'white', colour = 'white'),
  #      panel.background = element_rect(fill = 'white', colour = 'white'),
  #      axis.title = element_blank()
  )

library(patchwork)

micro_map + hylid_map + myo_map

# R script to generate the map wih the locations and the samples

```R
library(tidyverse)
library(dplyr)
require(maps)
library(ggrepel)
library(patchwork)
library(ggsci)
library(ggpubr)

# load world data
world_map <- map_data("world")

# load metadata
data <- read.delim("region_data.txt", sep="\t", header=T) %>%
  select("N", "country", "region", "region_code", "LAT", "LONG")

data$region <- str_replace(data$region, 'Sydney', 'New South Wales')
data$region_code <- str_replace(data$region_code, 'SYD', 'NSW')
data$LAT <- as.numeric(data$LAT)
data$LONG <- as.numeric(data$LONG)


# And now select data for USA map
states_map <- map_data("state")
data2 <- data %>% filter(country == 'USA')

#Preparing data for plotting

#Generating levels for population and super-population variables
data$region_code <- factor(data$region_code, 
                          levels = c('QUE', 'NSW', 
                                     'PAV',
                                     'SIC',
                                     'MCH', 'ILL', 'TEN', 'ARK', 
                                     'GEO', 'MIP', 'TEX', 'LOU'))


scale_colour_javier_PCA <- function(...){
  ggplot2:::manual_scale(
    'colour', 
    values = setNames(
      c('royalblue1', 'blue3',
        'turquoise3',
        'green',
        'lightgoldenrod1', 'darkgoldenrod1', 'orange', 'orange2',
        'orange3', 'tomato1', 'red2', 'darkred'), 
      c('QUE', 'NSW', 
        'PAV',
        'SIC',
        'MCH', 'ILL', 'TEN', 'ARK', 
        'GEO', 'MIP', 'TEX', 'LOU')), 
    ...
  )
}

# Make the maps
map <- ggplot() +
  geom_polygon(data = world_map, aes(x = world_map$long, y = world_map$lat, group = world_map$group), fill="grey90") +
  geom_point(data = data, aes(x = LONG, y = LAT, colour = region_code), size=2) +
  geom_text_repel(data = data, aes(x = LONG, y = LAT, label = paste0(region," (n = ",N,")")), size=4 , max.overlaps = Inf) +        
  theme_void() + 
  theme(legend.position = 'none')+
  ylim(-55, 90) +
  labs(colour="", shape="")+
  scale_colour_javier_PCA()

# save it
ggsave("Fig2a_worldmap_samplingsites.png", height=4, width=9)
```

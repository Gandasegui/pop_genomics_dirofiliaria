# R script to generate the map wih the locations and the samples

```R
library(ggplot2)
library(dplyr)
require(maps)
library(ggrepel)
library(patchwork)
library(ggsci)
library(ggpubr)

#
setwd("~/R/diro_newgenome/worldmap")
# load world data
world_map <- map_data("world")

# load metadata
data <- read.delim("region_data.txt", sep="\t", header=T) %>%
  select("N", "country", "region", "region_code", "LAT", "LONG")

data$region <- str_replace(data$region, 'Sydney', 'New South Wales')
data$region_code <- str_replace(data$region_code, 'SYD', 'NSW')
data[nrow(data) + 1,] <- c(1, 'China', 'Sichuan', "SIC", 30.813051538095294, 103.22450766942613)
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
a <- ggplot() +
  geom_polygon(data = world_map, aes(x = world_map$long, y = world_map$lat, group = world_map$group), fill="grey90") +
  geom_point(data = data, aes(x = LONG, y = LAT, colour = region_code), size=3) +
  geom_text_repel(data = data, aes(x = LONG, y = LAT, label = paste0(region," (n = ",N,")")), size=3.5) +        
  theme_void() + 
  labs(colour="", shape="") +
  ylim(-55,85) +
  xlim(-21,175) + 
  scale_colour_javier_PCA()


b <- ggplot() +
  geom_polygon(data = states_map, aes(x = states_map$long, y = states_map$lat, group = states_map$group), fill="grey90") +
  geom_point(data = data2, aes(x = LONG, y = LAT, colour = region_code), size=3) +
  geom_text_repel(data = data2, aes(x = LONG, y = LAT, label = paste0(region," (n = ",N,")")), size=3.5) +        
  theme_void() +
  labs(colour="", shape="") + 
  scale_colour_javier_PCA()

# And arrange them
ggarrange(b, a, labels = c("a", "b"), legend = 'none', hjust = -1, vjust = 1)

# save it
ggsave("worldmap_samplingsites.png", height=4, width=11)
```

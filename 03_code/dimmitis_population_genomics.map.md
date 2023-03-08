# Dirofilaria immitis populaiton genomics: sample sites

### stephen doyle




```bash
cd /nfs/users/nfs_s/sd21/lustre_link/dirofilaria_immitis/POPGEN/MAP
```


```R
library(tidyverse)
require(maps)
library(ggrepel)
library(patchwork)
library(ggsci)


world_map <- map_data("world")

# load metadata
data <- read.delim("map_metadata.txt", sep="\t", header=T)


region_colours <-
     c("PAV" = "#03C5CD",
     "NSW" = "#0000CD",
     "QLD" = "#4976FF",
     "TEN" = "#FFA503",
     "TEX" = "#EF0000",
     "GEO" = "#CD8500",
     "ARK" = "#EE9A00",
     "MCH" = "#FFEC8B",
     "LOU" = "#8B0000",
     "MIP" = "#FF6447",
     "ILL" = "#FFB90D")


# make a map
ggplot() +
     geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill="grey90") +
     geom_point(data = data, aes(x = longitude, y = latitude, colour = region_code, shape = country_code), size=3) +
     geom_text_repel(data = data, aes(x = longitude, y = latitude, label = paste0(region," (n = ", sample_n,")")), size=3, max.overlaps = Inf) +
     theme_void() +
     labs(title="A", colour="Region", shape="Country") +
     scale_colour_manual(values = region_colours, breaks=c('NSW', 'QLD', 'PAV', 'ARK', 'GEO', 'ILL', 'LOU', 'MCH', 'MIP', 'TEN', 'TEX'))

ggsave("diro_samples_map.pdf", height=5, width=7, useDingbats=FALSE)
ggsave("diro_samples_map.png")

```

![](../04_analysis/diro_samples_map.png)
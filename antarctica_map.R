# script for making map

# setwd("Google Drive/Projects/phaeo-miao/")

#packages to load
library(maps)
library(mapdata)
library(bitops)
library(RCurl)
library(png)
library(RJSONIO)
library(RgoogleMaps)
library(TeachingDemos)
library(dplyr)
library(mapproj)
library(rgdal)
library(maptools)
library(ggplot2)

## select out just the continent
antarctica <-  wrld_simpl[wrld_simpl$NAME == "Antarctica", ]

## define a sensible projection
# pr <- "+proj=laea +lat_0=-90 +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
pr <- CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
antarctica.laea <- spTransform(antarctica, CRS(pr))

ant.points <- fortify(antarctica.laea)

ant.points2 <- ant.points[which(ant.points$long != 0),]

ggplot(ant.points, aes(long, lat)) + geom_polygon()

ant_map <- ggplot(ant.points2, aes(long, lat, group = group)) + 
  geom_polygon(colour = 'grey40', fill = "grey60") + 
  theme(panel.background = element_rect(fill = "white"),
        plot.margin = margin(2, 2, 2, 2, "cm"), 
        plot.background = element_rect(fill = "white", colour = "black", size = 1),
        axis.ticks = element_blank(),
        axis.line = element_line(),
        axis.text = element_blank(),
        axis.title = element_blank()) + panel_border("black") + 
  annotate(geom = "point", x = 350000, y = -1340000, size = 8, colour = 'grey15', alpha = 0.9)

ggsave(ant_map, filename = "antarctica_map_mcmurdo.jpg")




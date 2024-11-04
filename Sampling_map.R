#######################################################################################
############################### Generate the sample map ##############################
#######################################################################################

# set the working directory
setwd("/Dinopium_GBS_on_laptop/Samping maps/")

# load the packages
library(maps)
library(mapdata)
library(maptools)
library(scales)
library(dplyr)
library(sf)

# read the coordinate data file
samps <- read.csv("sampling locations_2021 and 2018_decimal_coordinates.data.csv", stringsAsFactors = T)

## sub setting data according to the plumage color
samps_red <- filter(samps, pop_group_III == "Dino_R")
samps_hybrid <- filter(samps, pop_group_III == "Dino_H")
samps_yellow <- filter(samps, pop_group_III == "Dino_Y" | pop_group_III == "Dino_Y_J" | pop_group_III ==  "Dino_Y_M")


## Red the shp range data
Mannar <- st_read("mygeodata_Mannar/Mannar_peninsular-polygon.shp")
Jaffna <- st_read("Jaffna_Polygon/Jaffna.shp")

## Make the plot 
quartz(width = 7, height = 7)
map('worldHires', xlim = c(78, 83), ylim = c(5, 11), fill = T, col = "gray90")
plot(Mannar, add=T, col=alpha("blue", 0.6), border = F)
plot(Jaffna, add=T, col=alpha("green", 0.6), border = F)
points(79.82411,10.2845, pch=21, col="black", bg="gray3", cex=0.9)
text(80.4, 10.3,  "Point Calimere", cex=0.7)
points(samps_red$Long.decimal, samps_red$Lat.decimal, pch=21, col="black", bg="red", cex = 0.9) 
points(samps_hybrid$Long.decimal, samps_hybrid$Lat.decimal, pch=21, col="black", bg="orange", cex = 0.9) 
points(samps_yellow$Long.decimal, samps_yellow$Lat.decimal, pch=21, col="black", bg="yellow", cex = 0.9) 
map.axes(cex.axis=0.7)
map.scale(81.4, 5.3, ratio = F, relwidth = 0.2, cex=0.6)
title("The distribution of sampling loctions (n=140)", cex=1)


## Make the plot with jitter to the points on the map to avoid overlap.
# The jitter() function adds a small amount of random variation to each point's longitude and latitude, helping to visually separate overlapping points.
quartz(width = 7, height = 7)
map('worldHires', xlim = c(78, 83), ylim = c(5, 11), fill = T, col = "gray90")
plot(Mannar, add=T, col=alpha("blue", 0.6), border = F)
plot(Jaffna, add=T, col=alpha("green", 0.6), border = F)
points(79.82411,10.2845, pch=21, col="black", bg="gray3", cex=0.9)
text(80.4, 10.3,  "Point Calimere", cex=0.7)
points(jitter(samps_red$Long.decimal, amount = 0.05), 
       jitter(samps_red$Lat.decimal, amount = 0.05), 
       pch = 21, col = "black", bg = "red", cex = 0.9)
points(jitter(samps_hybrid$Long.decimal, amount = 0.05), 
       jitter(samps_hybrid$Lat.decimal, amount = 0.05), 
       pch = 21, col = "black", bg = "orange", cex = 0.9)
points(jitter(samps_yellow$Long.decimal, amount = 0.05), 
       jitter(samps_yellow$Lat.decimal, amount = 0.05), 
       pch = 21, col = "black", bg = "yellow", cex = 0.9)
map.axes(cex.axis=0.7)
map.scale(81.4, 5.3, ratio = F, relwidth = 0.2, cex=0.6)
title("The distribution of sampling loctions (n=140)", cex=1)


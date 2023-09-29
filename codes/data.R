rm(list = ls())

#import map
library(sf)
Rwanda_map <- st_read(".\\data\\RWA_adm\\RWA_adm0.shp")
plot(st_geometry(Rwanda_map))
st_crs(Rwanda_map)

Rwanda_map.web <-st_transform(Rwanda_map,3857)


#import data
data <- read.csv(".\\data\\Rwanda2020_individual_data.csv")
head(data)

ls()


#saveRDS(ASCARIS.estim,".\\outputs\\ASCARIS.2008.estim.rds")
#saveRDS(S.samples,".\\outputs\\ASCARIS.2008.S.samples.rds")




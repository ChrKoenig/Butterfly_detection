library(tidyverse)
library(png)

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")
rm(list = ls())

# ------------------------------------------------------------------------------------- #
# This script was developed by Christian König and has been tested on the four sample images in /Data/raw_data/traits/Zeuss. These images are named following the scheme 
# <Family><Genus><species><side of the photographed individual><sex of the photographed individual>.png
# The extraction of color characteristics for all species was done by Dirk Zeuss at Marburg University, who was hesitant to share the full dataset of scans.
# The results files can be found in the same folder (colors_extr_<family_name>.Rdata)

#### Process images ####
lookup = read_csv("~/Documents/lookup.csv")
pattern = paste(str_replace(lookup$name_Zeuss, " ", "_"), collapse = "|")
file_names = dir("Data/raw_data/traits/Zeuss/", recursive = T, full.names = T, pattern = pattern)

colors_extr = bind_rows(lapply(file_names, function(file_name){
  pic = readPNG(file_name, native = F)
  
  # Get mask
  alpha = as.vector(pic[,,4]) 
  
  # extract color infos
  red = as.vector(pic[,,1])[alpha != 0] 
  green = as.vector(pic[,,2])[alpha != 0]
  blue = as.vector(pic[,,3])[alpha != 0]
 
  # convert to HSV (hue, saturation, value) color space 
  pic_hsv = rgb2hsv(red, green, blue, maxColorValue = 1)
  mean_sat = mean(pic_hsv[2,])
  mean_lgt = mean(pic_hsv[3,])
  if(mean_sat < 0.3){ 
    main_color = "none"   # no apparent color
  } else {
    hue_segmented = cut(pic_hsv[1,]*360, breaks = c(0, 20, 50, 70, 160, 200, 280, 330, 361), include.lowest=TRUE, right=FALSE,
                             labels = c("red", "orange", "yellow", "green", "cyan", "blue", "magenta", "red"))
    color_freq = table(factor(hue_segmented))
    main_color = names(color_freq[which.max(color_freq)]) # most frequent color
  }
  
  file_name_short = gsub("^(.*/)(.*)(\\.png$)", "\\2", file_name)
  genus = strsplit(file_name_short, split = "_")[[1]][2]
  epithet = strsplit(file_name_short, split = "_")[[1]][3]
  species = paste(genus, epithet)
  side = ifelse(grepl("_t_", file_name_short), "top", "bottom")
  sex = ifelse(grepl("_m$", file_name_short), "male", "female")
  
  return(tibble(species = species, side = side, sex = sex, main_color = main_color, mean_sat = mean_sat, mean_lgt = mean_lgt))
}))
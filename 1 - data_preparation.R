library("raster")
library("sp")
library("tidyverse")
library("taxize")

setwd("~/ownCloud/Projects/Berlin/06 - Butterfly_detection/")
rm(list = ls())

# ------------------------------------------------------------------------------------- #
#### Species data ####
butterflies = read_csv("Data/raw_data/1550_TF_Daten_blinded.csv",
                       col_types = list(ArtCodeTagf = "f", NameTagfLat = "f", aldStao_ = "f", Jahr_ = "f", Bearbeiter = "f")) %>% 
  filter(Aufnahmetyp == "Normalaufnahme_Z7", Flags.ZA == 0) %>% 
  na_if("?") %>% 
  mutate_at(vars(PreAbs1:PreAbs7), function(x){
    x = replace(x, x == -1, NA) # Not visited
    x = replace(x, x == 0, 0) # Not observed 
    x = replace(x, x == 1 | x == 2, 1) # Observed once
    x = replace(x, x > 2, 2)  # Observed twice
  }) %>% 
  mutate_at(vars(contains("Day")), list(~ as.numeric(.))) %>% 
  mutate(region = case_when(BGR %in% c("Östliche Zentralalpen","Westliche Zentralalpen") ~ "ZA",
                            BGR %in% c("Jura") ~ "JU",
                            BGR %in% c("Mittelland") ~ "ML",
                            BGR %in% c("Alpennordflanke") ~ "AN",
                            BGR %in% c("Alpensüdflanke") ~ "AS"),
         species = as.factor(str_replace(NameTagfLat, "-Komplex", "")),
         coord_x = coordx * 1000, 
         coord_y = coordy * 1000)  %>% 
  rename_all(tolower) %>% 
  dplyr::select(spec_id = artcodetagf, species = nametagflat, site_id = aldstao_, region, coord_x, coord_y, 
                elevation = hoehe, year = jahr_, observer = bearbeiter, day1:day7, preabs1:preabs7)

save(butterflies, file = "Data/butterflies.RData") 

# ------------------------------------------------------------------------------------- #
#### Sampling sites ####
site_coords = butterflies %>%  
  dplyr::select(site_id, coord_x, coord_y) %>% 
  distinct()

site_polygons = lapply(site_coords$site_id, function(i){
  x_i = site_coords$coord_x[site_coords$site_id == i]
  y_i = site_coords$coord_y[site_coords$site_id == i]
  coords_i = matrix(c(x_i - 500, y_i - 500, 
                      x_i - 500, y_i + 500,
                      x_i + 500, y_i + 500, 
                      x_i + 500, y_i - 500, 
                      x_i - 500, y_i - 500), ncol = 2, byrow = T)
  p_i = Polygon(coords_i)
  ps_i = Polygons(list(p_i), ID = i)
  sps_i = SpatialPolygons(list(ps_i))
  return(sps_i)
})

site_polygons = do.call("rbind", c(args = site_polygons, makeUniqueIDs = TRUE))
extent(site_polygons)
plot(site_polygons)

# ------------------------------------------------------------------------------------- #
#### Env data ####
climtopo = stack("Data/raw_data/CH_climtopo_1km.grd")
landuse = stack("Data/raw_data/Areal_NOAS04_prop1km.grd")
landuse_key = read_csv("Data/raw_data/landuse_types.csv")

extent(climtopo)
extent(landuse)

# Extract values
extr_climtopo = raster::extract(climtopo, site_polygons, fun = median)
extr_landuse  = raster::extract(subset(landuse, grep("AS09_17", names(landuse))), site_polygons, fun = median) # TODO: Discuss appropriate level with damaris, Some rows don't add up tp 100

# Create SPDF
extr_env = cbind(extr_climtopo, extr_landuse) %>% 
  as_tibble() %>% 
  mutate(site_id = as.integer(names(site_polygons))) %>% 
  left_join(distinct(butterflies[,c("site_id", "elevation")])) %>% 
  column_to_rownames("site_id")

sample_sites = SpatialPolygonsDataFrame(site_polygons, extr_env, match.ID = T)
save(sample_sites, file =  "Data/sample_sites.RData")

# ------------------------------------------------------------------------------------- #
#### Trait data ####
traits_altermatt = read_csv("Data/raw_data/traits/traits_Altermatt.csv") %>% 
  select(species, generations, hibernation, larva_diet)

traits_zeuss = read_csv("Data/raw_data/traits/bodysize_and_lightness_Europe.csv") %>% 
  rowwise() %>% 
  mutate(genus = strsplit(file, split = "_")[[1]][2],
         epithet = strsplit(file, split = "_")[[1]][3],
         face = ifelse(grepl("_t_", file), "top", "bottom"),
         sex = ifelse(grepl("_m.png", file), "male", "female"),
         species = paste(genus, epithet)) %>% 
  dplyr::select(species, meanRGB , body_area = bodysize_cm2, face, sex) %>% 
  pivot_wider(names_from = face, names_prefix = "RGB_", values_from = meanRGB) %>% 
  group_by(species) %>% 
  dplyr::summarise(body_area = mean(body_area, na.rm = T), 
                   RGB_bottom = mean(RGB_bottom, na.rm = T),
                   RGB_top = mean(RGB_top, na.rm = T))

traits_final = full_join(traits_altermatt, traits_zeuss)
save(traits_final, file =  "Data/traits_final.RData")

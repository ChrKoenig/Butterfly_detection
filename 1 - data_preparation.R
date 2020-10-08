library("raster")
library("sp")
library("tidyverse")

setwd("~/Butterfly_project/")
rm(list = ls())

# ------------------------------------------------------------------------------------- #
#### Species data ####
observations = read_csv("Data/raw_data/1550_TF_Daten_blinded.csv",
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
  dplyr::select(spec_id = artcodetagf, species, site_id = aldstao_, region, coord_x, coord_y, 
                elevation = hoehe, year = jahr_, observer = bearbeiter, day1:day7, preabs1:preabs7)

save(observations, file = "Data/observations.RData") 

# obtain unique combinations of site_id x observer x year
visit_info = observations %>% 
  dplyr::select(site_id:observer, day1:day7) %>% 
  distinct()

# Complete observation data, i.e. add absences for unobserved species 
observations_completed = observations %>% 
  filter(as.numeric(as.character(year)) > 2007) %>% # No good observations before 2007
  droplevels() %>% 
  dplyr::select(spec_id, site_id, year, preabs1:preabs7) %>% 
  complete(spec_id, site_id, year) %>% 
  replace(is.na(.), 0) %>% 
  inner_join(visit_info, by = c("site_id", "year")) %>% 
  inner_join(distinct(select(observations, "species", "spec_id"))) %>% 
  relocate(species, .after = spec_id) %>% 
  relocate(preabs1:preabs7, .after = last_col()) %>% 
  mutate_if(is.factor, as.character) 

save(observations_completed, file = "Data/observations_completed.RData")
# ------------------------------------------------------------------------------------- #
#### Sampling sites ####
site_coords = observations %>%  
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
  mutate(site_id = factor(names(site_polygons))) %>% 
  left_join(distinct(observations[,c("site_id", "elevation")])) %>% 
  column_to_rownames("site_id")

sample_sites = SpatialPolygonsDataFrame(site_polygons, extr_env, match.ID = T)
save(sample_sites, file =  "Data/sample_sites.RData")

# ------------------------------------------------------------------------------------- #
#### Trait data ####
# traits_altermatt = read_csv("Data/raw_data/traits/traits_Altermatt.csv") %>% 
#  select(species, generations, hibernation, larva_diet)

synonym_lookup = drop_na(read_csv("Data/synonym_lookup_edit.csv"))

traits_zeuss = read_delim("Data/raw_data/traits/Europa_schmetterlinge_data_art_final.csv", delim = ";", locale = locale(decimal_mark = ",")) %>% 
  rowwise() %>% 
  mutate(genus = strsplit(file, split = "_")[[1]][2],
         epithet = strsplit(file, split = "_")[[1]][3],
         species = paste(genus, epithet),
         color = case_when(meanHue*360 < 20 ~ "red",
                           meanHue*360 < 50 ~ "orange",
                           meanHue*360 < 70 ~ "yellow",
                           meanHue*360 < 160 ~ "green",
                           meanHue*360 < 200 ~ "cyan",
                           meanHue*360 < 280 ~ "blue",
                           meanHue*360 < 330 ~ "magenta",
                           TRUE ~ "red")) %>% 
  full_join(synonym_lookup, by = c("species" = "name_Zeuss")) %>% 
  mutate(species = replace(species, !is.na(name_BDM), name_BDM)) %>% 
  dplyr::select(species, generations = voltinism_num, color, saturation = meanChroma, 
                brightness = meanValue, body_area = bodysize_cm2) 

traits_final = traits_zeuss
save(traits_final, file =  "Data/traits_final.RData")

# ------------------------------------------------------------------------------------- #
#### Final species selection (observed + traits available)
# traits_only = sort(setdiff(traits_final$species, unique(observations_completed$species)))
# obs_only = sort(setdiff(unique(observations_completed$species), traits_final$species))
# 
# # Manually match synonyms, BDM naming convention has priority
# write_csv(data.frame(name_BDM = obs_only, name_Zeuss = ""), path = "Data/synonym_lookup.csv")
# 
# # Load species lookup table
# synonym_lookup = drop_na(read_csv("Data/synonym_lookup_edit.csv"))

species_final = traits_final %>% 
  inner_join(distinct(select(observations_completed, "spec_id", "species"))) %>% 
  select(species, spec_id)

save(species_final, file = "Data/species_final.RData")

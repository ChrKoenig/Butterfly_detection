library("raster")
library("sp")
library("tidyverse")

setwd("~/Butterfly_project/")
rm(list = ls())

# ------------------------------------------------------------------------------------- #
#### Species data ####
observations = read_csv("Data/raw_data/1550_TF_Daten_blinded.csv",
                       col_types = list(ArtCodeTagf = "f", NameTagfLat = "f", aldStao_ = "f", Jahr_ = "n", Bearbeiter = "f")) %>% 
  filter(Aufnahmetyp == "Normalaufnahme_Z7", Flags.ZA == 0, Jahr_ >= 2007) %>% 
  na_if("?") %>% 
  mutate_at(vars(PreAbs1:PreAbs7), function(x){
    x = replace(x, x == -1, NA) # Not visited
    x = replace(x, x == 0, 0) # Not observed 
    x = replace(x, x == 1 | x == 2, 1) # Observed once
    x = replace(x, x > 2, 2)  # Observed twice
  }) %>% 
  mutate_at(vars(contains("Day")), list(~ as.numeric(.))) %>% 
  mutate(species = as.factor(str_replace(NameTagfLat, "-Komplex", "")),
         coord_x = coordx * 1000, 
         coord_y = coordy * 1000,
         year = as.factor(Jahr_))  %>% 
  rename_all(tolower) %>% 
  dplyr::select(spec_id = artcodetagf, species, site_id = aldstao_, coord_x, coord_y, 
                elevation = hoehe, year = jahr_, observer = bearbeiter, day1:day7, preabs1:preabs7)

save(observations, file = "Data/observations.RData") 

# obtain unique combinations of site_id x observer x year
visit_info = observations %>% 
  dplyr::select(site_id:observer, day1:day7) %>% 
  distinct()

# Count observations per species / year*site (ignore repeated sampling)
spec_count = table(observations[rowSums(observations[,paste0("preabs", 1:7)], na.rm = T) > 0, "spec_id"])
observations_count = spec_count %>% 
  as.data.frame() %>% 
  rename("spec_id" = "Var1", "n" = "Freq") %>% 
  filter(n > 30) 

# Complete observation data, i.e. add absences for unobserved species 
observations_completed = observations %>% 
  filter(spec_id %in% observations_count$spec_id) %>% # Use only species with > 30 observations
  droplevels() %>% 
  dplyr::select(spec_id, site_id, year, preabs1:preabs7) %>% 
  complete(spec_id, site_id, year) %>% 
  replace(is.na(.), 0) %>% 
  inner_join(visit_info, by = c("site_id", "year")) %>% 
  inner_join(distinct(dplyr::select(observations, "species", "spec_id"))) %>% 
  relocate(species, .after = spec_id) %>% 
  relocate(preabs1:preabs7, .after = last_col()) %>% 
  mutate_if(is.factor, as.character) 

save(observations_completed, file = "Data/observations_completed.RData")

# ------------------------------------------------------------------------------------- #
#### Final species selection (observed + traits available)
# sort(setdiff(unique(observations_completed$species), colors_extr$species))
# write_csv(data.frame(name_BDM = sort(setdiff(unique(observations_completed$species), colors_extr$species)), name_Zeuss = NA),
#          "Data/raw_data/traits/Zeuss/synonym_lookup.csv")
# 
# # Load species lookup table
synonyms_zeuss = read_csv("Data/raw_data/traits/Zeuss/synonym_lookup_edit.csv")
synonyms_embt = read_csv("Data/raw_data/traits/EMBTv1.2/synonym_lookup_edit.csv")
species_final = distinct(dplyr::select(observations_completed, "spec_id", "species")) %>% 
  anti_join(filter(synonyms_zeuss, is.na(name_Zeuss)), by = c("species" = "name_BDM")) %>% 
  anti_join(filter(synonyms_embt, is.na(name_EMBT)), by = c("species" = "name_BDM")) %>% 
  dplyr::select(species, spec_id)

save(species_final, file = "Data/species_final.RData")

# ------------------------------------------------------------------------------------- #
#### Trait data ####
# colors_extr = bind_rows(lapply(list.files(path = "Data/raw_data/traits/Zeuss/", pattern = "idae\\.Rdata", full.names = T), function(file){
# +     load(file)
# +     return(colors_extr)
# + }))
# save(colors_extr, file = "Data/raw_data/traits/Zeuss/colors_extr.RData")
# Color traits: main color, saturation, lightness
synonyms_zeuss = drop_na(read_csv("Data/raw_data/traits/Zeuss/synonym_lookup_edit.csv"))
load("Data/raw_data/traits/Zeuss/colors_extr.RData")
traits_zeuss = colors_extr %>% 
  group_by(species, side) %>% 
  summarize(main_color = raster::modal(main_color),
            mean_sat = mean(mean_sat),
            mean_lgt = mean(mean_lgt)) %>% 
  pivot_wider(id_cols = species, names_from = side, values_from = c(main_color, mean_sat, mean_lgt)) %>% 
  left_join(synonyms_zeuss, by = c("species" = "name_Zeuss")) %>%
  rowwise() %>% 
  mutate(species = replace(species, !is.na(name_BDM), name_BDM)) %>% 
  filter(species %in% species_final$species) %>% 
  dplyr::select(-name_BDM)

# Morphological/behavioural traits: Voltinism, Wing index, Hostplant specificity index, Basking Site
synonyms_embt = drop_na(read_csv("Data/raw_data/traits/EMBTv1.2/synonym_lookup_edit.csv"))
traits_embt = read_csv("Data/raw_data/traits/EMBTv1.2/EMBT_trait_states.csv") %>% 
  mutate(species = str_replace(Taxon, "_", " ")) %>% 
  dplyr::select("species", "Vol_min", "WIn", "HSI", "FMo_Average") %>% 
  full_join(synonyms_embt, by = c("species" = "name_EMBT")) %>% 
  rowwise() %>% 
  mutate(species = replace(species, !is.na(name_BDM), name_BDM)) %>% 
  dplyr::select(-name_BDM)

# Join traits
traits_final = traits_embt %>% 
  full_join(traits_zeuss) %>% 
  dplyr::filter(species %in% species_final$species)

save(traits_final, file =  "Data/traits_final.RData")

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
  inner_join(distinct(dplyr::select(observations_completed, "spec_id", "species"))) %>% 
  dplyr::select(species, spec_id)

save(species_final, file = "Data/species_final.RData")

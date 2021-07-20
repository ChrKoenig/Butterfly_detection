library("runjags") 

# ------------------------------------------------------- #
#              Sample DIC for fitted models            ####
# ------------------------------------------------------- #
MSOM1 = readRDS("//import/calc9z/data-zurell/koenig/Butterfly_models_fit/MSOM1/MSOM1_15.RDS")              # Read last extension of fitted models
MSOM1_dic = runjags::extract(MSOM1, "DIC")                                                                 # Extract DIC
save(MSOM1_dic, file = "//import/calc9z/data-zurell/koenig/Butterfly_models_fit/MSOM1/MSOM1_dic.RData")    # save result

MSOM2 = readRDS("//import/calc9z/data-zurell/koenig/Butterfly_models_fit/MSOM2/MSOM2_15.RDS")              # Read last extension of fitted models     
MSOM2_dic = runjags::extract(MSOM2, "DIC")                                                                 # Extract DIC  
save(MSOM2_dic, file = "//import/calc9z/data-zurell/koenig/Butterfly_models_fit/MSOM2/MSOM2_dic.RData")    # save result

MSOM3 = readRDS("//import/calc9z/data-zurell/koenig/Butterfly_models_fit/MSOM3/MSOM3_15.RDS")              # Read last extension of fitted models 
MSOM3_dic = runjags::extract(MSOM3, "DIC")                                                                 # Extract DIC    
save(MSOM3_dic, file = "//import/calc9z/data-zurell/koenig/Butterfly_models_fit/MSOM3/MSOM3_dic.RData")    # save result 


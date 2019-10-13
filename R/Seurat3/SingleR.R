library(SingleR)
source("../R/SingleR_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/Lymphoma.data_EC_10_20190922.Rda"))
ref.se <- BlueprintEncodeData(rm.NA = "rows")

singler = SingleR(test = object_data, ref = ref.se, 
                  labels = ref.se$label.main)

save(singler,file="output/singler_Lymphoma_EC_10T_20190922.Rda")

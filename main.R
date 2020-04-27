#' Automatic Topographic Correction  ---------------------------------------
#' @author Cesar and Lesly
#' @Methodology Sun Canopy Sensor + C correction
#' @paper https://www.mdpi.com/2220-9964/6/9/287/pdf

library(rgee)
library(sf)
source("functions.R")
options(rgee.print.option = "simply")

ee_reattach()
ee_Initialize(drive = TRUE)

# Region of Interest
rajucolta <- st_read("data/rajucolta.geojson") %>% sf_as_ee()

# 1. Move From Earth Engine dataset to EE asset
s2_raw <- get_s2_img(c("2019-09-16", "2019-09-17"))

# 2. Remove Clouds
s2_nocloud <- maskS2clouds(s2_raw)

# 3. Topographic Correction (SCS + C)
s2_topo <- topographic_correction_s2(s2_img = s2_nocloud)

# 4. Download Results
to_download <- ee$ImageCollection(c(s2_raw$toInt(), s2_nocloud$toInt(), s2_topo$toInt()))
ee_imagecollection_to_local(
  ic = to_download,
  region = rajucolta$geometry(),
  dsn = paste0("output/",c("raw.tif","nocloud.tif","topo.tif")),
  scale = 10,
  via = "drive"
)

# 5. Display results in mapview
Map$centerObject(s2_raw)
Map$addLayer(s2_raw, list(min = 0, max = 10000, bands = c("B3","B2","B1")), "s2_raw") +
Map$addLayer(s2_nocloud, list(min = 0, max = 10000, bands = c("B3","B2","B1")), "s2_nocloud") +
Map$addLayer(s2_topo, list(min = 0, max = 10000, bands = c("2_RED","1_GREEN","0_BLUE")), "s2_ic")

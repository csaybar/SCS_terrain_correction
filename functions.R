# Get S2 image
get_s2_img <- function(date = c("2019-09-16", "2019-09-17")) {
  s2_raw <- ee$ImageCollection("COPERNICUS/S2_SR")$
    filterBounds(rajucolta)$
    filterDate(date[1], date[2])$
    first()
  image_id <- paste0(ee_get_assethome(), "/", "s2_rajucolta")
  task <- ee_image_to_asset(
    image = s2_raw,
    overwrite = TRUE,
    region = rajucolta$geometry(),
    assetId = image_id
  )
  task$start()
  ee_monitoring(task)
  ee$Image(image_id)
}

# Function to mask clouds using the Sentinel-2 QA band.
maskS2clouds <- function(image, range = c(10, 10)) {
  qa <- image$select("SCL")$divide(0.0001) # warning this mask could change
  mask <- qa$gte(range[1])$And(qa$lte(range[2]))$eq(0)

  # Return the masked and scaled data, without the QA bands.
  img_clean <- ee$Image(
    image$updateMask(mask)$
      select("B.*")$
      copyProperties(image, list("system:time_start"))
  )
  image_id <- paste0(ee_get_assethome(), "/", "s2_clean_rajucolta")
  task <- ee_image_to_asset(
    image = img_clean,
    overwrite = TRUE,
    assetId = image_id,
    region = rajucolta$geometry()
  )
  task$start()
  ee_monitoring(task)
  ee$Image(image_id)
}

#' SCS terrain correction
apply_SCSccorr <- function(band) {
  mask1 <- img$select("NIR")$gt(-0.1)
  mask2 <- img$select("slope")$gte(5)$
    And(img$select("IC")$gte(0))$
    And(img$select("NIR")$gt(-0.1))
  img_2 <- ee$Image(img$updateMask(mask2))

  out <- img_2$select("IC", ee$String(band))$reduceRegion(
    reducer = ee$Reducer$linearFit(), # Compute coefficients: a(slope), b(offset), c(b/a)
    geometry = ee$Geometry(img$geometry()$buffer(-5000)), # trim off the outer edges of the image for linear relationship
    scale = 10,
    maxPixels = 1000000000
  )
  out_a <- ee$Number(out$get("scale"))
  out_b <- ee$Number(out$get("offset"))
  out_c <- ee$Number(out$get("offset"))$divide(ee$Number(out$get("scale")))

  # apply the SCSc correction
  SCSc_output <- img_2$expression("((image * (cosB * cosZ + cvalue)) / (ic + cvalue))", list(
    "image" = img_2$select(ee$String(band)),
    "ic" = img_2$select("IC"),
    "cosB" = img_2$select("cosS"),
    "cosZ" = img_2$select("cosZ"),
    "cvalue" = out_c
  ))
  ee$Image(SCSc_output)
}

topoCorr_IC_L4toL8 <- function(img) {
  dem <- ee$Image("USGS/SRTMGL1_003")
  # Extract image metadata about solar position
  SZ_rad <- ee$Image$constant(ee$Number(img$get("MEAN_SOLAR_ZENITH_ANGLE")))$multiply(pi)$divide(180)$clip(img$geometry()$buffer(10000))
  SA_rad <- ee$Image$constant(ee$Number(img$get("MEAN_SOLAR_AZIMUTH_ANGLE"))$multiply(pi)$divide(180))$clip(img$geometry()$buffer(10000))

  # Creat terrain layers
  slp <- ee$Terrain$slope(dem)$clip(img$geometry()$buffer(10000))
  slp_rad <- ee$Terrain$slope(dem)$multiply(pi)$divide(180)$clip(img$geometry()$buffer(10000))
  asp_rad <- ee$Terrain$aspect(dem)$multiply(pi)$divide(180)$clip(img$geometry()$buffer(10000))

  # Calculate the Illumination Condition (IC)
  # slope part of the illumination condition
  cosZ <- SZ_rad$cos()
  cosS <- slp_rad$cos()

  slope_illumination <- cosS$expression(
    expression = "cosZ * cosS",
    opt_map = list(cosZ = cosZ, cosS = cosS$select("slope"))
  )
  # aspect part of the illumination condition
  sinZ <- SZ_rad$sin()
  sinS <- slp_rad$sin()
  cosAziDiff <- (SA_rad$subtract(asp_rad))$cos()
  aspect_illumination <- sinZ$expression(
    expression = "sinZ * sinS * cosAziDiff",
    opt_map = list(sinZ = sinZ, sinS = sinS, cosAziDiff = cosAziDiff)
  )
  # full illumination condition (IC)
  ic <- slope_illumination$add(aspect_illumination)

  # Add IC to original image
  img_plus_ic <- ee$Image(img$addBands(ic$rename("IC"))$
    addBands(cosZ$rename("cosZ"))$
    addBands(cosS$rename("cosS"))$
    addBands(slp$rename("slope")))
  return(img_plus_ic)
}

topoCorr_SCSc_L4toL8 <- function(img) {
  img <<- img
  mask1 <- img$select("NIR")$gt(-0.1)
  mask2 <- img$select("slope")$gte(5)$
    And(img$select("IC")$gte(0))$
    And(img$select("NIR")$gt(-0.1))
  img_2 <<- ee$Image(img$updateMask(mask2))

  # Specify Bands to topographically correct
  bandList <- ee$List(list("BLUE", "GREEN", "RED", "NIR", "SWIR1", "SWIR2"))
  img_SCSccorr <- ee$ImageCollection(bandList$map(ee_pyfunc(apply_SCSccorr)))$
    toBands()$
    addBands(img$select("IC"))
  bandList_IC <- ee$List(c(bandList, "IC"))$flatten()
  img_SCSccorr$unmask(img$select(bandList_IC))$
    addBands(mask1$rename("initMask"))$
    addBands(mask2$rename("corrMask"))
}

topographic_correction_s2 <- function(s2_img) {
  s2_img <- s2_img$select(
    opt_selectors = c("B2", "B3", "B4", "B8", "B11", "B12"),
    opt_names = c("BLUE", "GREEN", "RED", "NIR", "SWIR1", "SWIR2")
  )

  # Map$addLayer(s2_img,list(min = 0, max = 3000, bands = "SWIR1,NIR,RED"), "original")
  # Map$addLayer(s2_img, list(min = 0, max = 3000, bands = "SWIR1,NIR,RED"), "false color")
  s2_img %>%
    topoCorr_IC_L4toL8() %>%
    topoCorr_SCSc_L4toL8() ->
    s2_scs

  s2_topo <- s2_scs$select(c("0_BLUE", "1_GREEN", "2_RED", "3_NIR", "4_SWIR1", "5_SWIR2"))$
    rename(c("BLUE", "GREEN", "RED", "NIR", "SWIR1", "SWIR2"))

  image_id <- paste0(ee_get_assethome(), "/", "s2_terrain_rajucolta")
  task <- ee_image_to_asset(
    image = s2_topo,
    overwrite = TRUE,
    assetId = image_id,
    region = rajucolta$geometry()
  )
  task$start()
  ee_monitoring(task)
  ee$Image(image_id)
}

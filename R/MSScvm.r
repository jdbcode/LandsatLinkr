#' Create a cloud and cloud shadow mask for Landsat MSS imagery
#' 
#' Takes in any numeric value and squares it.
#' @param file MSS image filename (full system path to MSS file) 
#' @param demfile DEM filename (full system path to spatially coincident DEM file)
#' @return A binary raster with the same dimensions as the MSS image where pixels with value 1 represent clear pixel and 0 as obsured by either cloud or cloud shadow
#' @import raster
#' @import SDMTools
#' @export


msscvm = function(file, demfile){
  
#    library(SDMTools)
#    library(raster)
#    library(MSScvm)
# #   
#    file = "E:/mss/wrs1/036032/images/1975/LM20360321975226_radiance.tif"
#    demfile = "K:/gis_data/dems/wrs_dem/wrs1_036032_60m_dem.tif"
  
  #get the metadata
  info = get_metadata(file)
  
  b1 = as.matrix(raster(file, 1))
  b2 = as.matrix(raster(file, 2))
  b4 = as.matrix(raster(file, 4))
  
  #apply topographic correction to band 4 for identifying cloud shadows
  k=0.55
  sunzen = info$sunzen*(pi/180)
  
  #crop the hillshade layer and convert to illumination
  dem = raster(demfile)
  ref = raster(file)
  dem_ex  = alignExtent(dem, ref, snap="near")
  extent(dem) = dem_ex
  dem = crop(dem,ref)

#dem = crop(raster(demfile), raster(file))
  slope = terrain(dem, opt="slope")
  aspect = terrain(dem, opt="aspect")
  ill = as.matrix(hillShade(slope, aspect, angle=info$sunelev, direction=info$sunaz, normalize=F))
  
  #apply the correction
  c = (cos(sunzen)/ill)^k
  b4topoc = round(b4*c)
  
  #find clouds
  cloud = (b1 - b2)/(b1 + b2)
  cloudtight = which(cloud > 0.0450 & b1 > 7500 | b1 > 14000)
  
  #find the shadows  
  b4topoc[cloudtight] = NA
  
  b4nocldmean = mean(b4topoc, na.rm=T)
  shadowthresh1 = round(0.40*b4nocldmean + 536.94)
  nocldorshdw = which(b4topoc > shadowthresh1) #find non-shadow pixels firest pass
  b4nocldorshdw = b4topoc[nocldorshdw]
  b4nocldorshdwmean = mean(b4nocldorshdw, na.rm=T)
  shadowthresh2 = round(0.47*b4nocldorshdwmean + 158.56)
  shadows = which(b4topoc <= shadowthresh2)
  
  #find the water
  slope = as.matrix(slope)  
  ndvi = (b4-b2)/(b4+b2)
  waterindex = which(ndvi < -0.2 & slope < (0*(pi/180)))    #slope = 0.0 degrees
  shoreindex = which(ndvi < -0.3 & slope < (0.5*(pi/180)))  #slope = 0.5 degrees
  
  b2[] = 0
  b2[c(waterindex,shoreindex)] = 1
  clumps = .Call("ccl", b2, PACKAGE = "SDMTools")  
  clumps = setValues(raster(file), clumps)
  #clumps = ConnCompLabel(b2)
  
  fre = freq(clumps)
  these = which(fre[,2] < 4)
  values = fre[these,1]
  m = match(as.matrix(clumps), values) 
  these = which(is.na(m) == F)
  #b2 = as.matrix(b2)
  b2[these] = 0
  b2 = setValues(raster(file),b2)
  
  b2 = focal(b2, w=matrix(1,5,5), fun=max) 
  b2 = as.matrix(b2)
  
  #create a blank raster and fill it in
  b1[] = 0
  b1[cloudtight] = 1  
  b1[shadows] = 1
  b1[b2 == 1] = 0 #set water back to clear
  #b1 = setValues(raster(file), b1)
  clumps = .Call("ccl", b1, PACKAGE = "SDMTools")  
  clumps = setValues(raster(file), clumps)
  #clumps = ConnCompLabel(b1)  
  
  fre = freq(clumps)
  these = which(fre[,2] < 10)
  values = fre[these,1]
  m = match(as.matrix(clumps), values) 
  these = which(is.na(m) == F)
  #b1 = as.matrix(b1)
  b1[these] = 0
  b1 = setValues(raster(file),b1)
  
  b1 = focal(b1, w=matrix(1,5,5), fun=max)
  b1[] = as.numeric(values(b1) == 0) 
  b1[is.na(raster(file))] = NA
  projection(b1) = set_projection(file)
  b1 = as(b1, "SpatialGridDataFrame")
  outfile = sub("radiance", "cloudmask", file)
  writeGDAL(b1, outfile, drivername = "GTiff", type = "Byte", mvFlag=255, options="INTERLEAVE=BAND")
}






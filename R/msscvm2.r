#' Create a cloud and cloud shadow mask for Landsat MSS imagery
#' 
#' Takes in any numeric value and squares it.
#' @param file character. MSS reflectance image filename (full system path to MSS file) 
#' @param demfile character. DEM filename (full system path to spatially coincident DEM file)
#' @return A binary raster with the same dimensions as the MSS image where pixels with value 1 represent clear pixel and 0 as obsured by either cloud or cloud shadow
#' @import raster
#' @import rgdal
#' @import SDMTools
#' @export


msscvm2 = function(file, demfile, test=T){
  #msscvm
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
  clouds = which(cloud > 0.0 & b1 > 1750 | b1 > 3900)
  
  #find the shadows  
  b4topoc[clouds] = NA
  
  b4nocldmean = mean(b4topoc, na.rm=T)
  shadowthresh1 = round(0.40*b4nocldmean + 247.97) #rad = 536.94
  nocldorshdw = which(b4topoc > shadowthresh1) #find non-shadow pixels first pass
  b4nocldorshdw = b4topoc[nocldorshdw]
  b4nocldorshdwmean = mean(b4nocldorshdw, na.rm=T)
  shadowthresh2 = round(0.47*b4nocldorshdwmean + 73.23) #refl = 73.23 rad = 158.56
  shadows = which(b4topoc <= shadowthresh2)
  
  #find the water
  slope = as.matrix(slope)  
  ndvi = (b4-b2)/(b4+b2)
  waterindex = which(ndvi < 0.0850 & slope < (0.5*(pi/180))) #| ndvi < -0.05 & & slope < (1.5*(pi/180))
  
  #create a set of blank rasters for filtering
  b1[] = 0
  water=shadow=cloud=b1 
  b1=0
  
  water[waterindex] = 1
  clumps = .Call("ccl", water, PACKAGE = "SDMTools")  
  clumps = setValues(ref, clumps)
  fre = freq(clumps)
  these = which(fre[,2] < 7)
  values = fre[these,1]
  m = match(as.matrix(clumps), values) 
  these = which(is.na(m) == F)
  water[these] = 0
  
  water = setValues(ref,water)
  water = focal(water, w=matrix(1,5,5), fun=max)  
  waterindex = which(as.matrix(water) == 1)
  
  shadow[shadows] = 1 
  shadow[waterindex] = 0
  cloud[clouds] = 1
  
  #filter the aggregated cloud and shadow 
  clumps = .Call("ccl", cloud, PACKAGE = "SDMTools")  
  clumps = setValues(ref, clumps)
  fre = freq(clumps)
  these = which(fre[,2] < 10) #10
  values = fre[these,1]
  m = match(as.matrix(clumps), values) 
  these = which(is.na(m) == F)
  cloud[these] = 0
   
  cloud = setValues(ref,cloud)
  cloud = focal(cloud, w=matrix(1,5,5), fun=max)
  
  reso = xres(ref)
  if(reso == 60){
    r = raster(ncol=31,nrow=31)
    ext = extent(0, 31, 0, 31)
    extent(r) = ext
    projection(r) = set_projection(file)
    r[] = NA
    r[16,16] = 1
    dist = gridDistance(r, 1)
    kernal = dist <= 16
  }
  if(reso == 30){
    r = raster(ncol=61,nrow=61)
    ext = extent(0, 61, 0, 61)
    extent(r) = ext
    projection(r) = set_projection(file)
    r[] = NA
    r[31,31] = 1
    dist = gridDistance(r, 1)
    kernal = dist <= 31
  }
  
  shadowproj = focal(cloud, w=as.matrix(kernal), fun=max)
  shiftstart = 1000/tan((pi/180)*info$sunelev)
  shiftend = 7000/tan((pi/180)*info$sunelev)
    
  shiftlen = seq(shiftstart,shiftend,900)
    
  shiftit = function(shadowproj,shiftlen,info,reso){
    if(info$sunaz > 90 & info$sunaz < 180){
      angle = info$sunaz-90  
      yshift = round((sin((pi/180)*angle) * shiftlen)/reso)*reso
      xshift = round((cos((pi/180)*angle) * shiftlen * -1)/reso)*reso
      shadowproj = shift(shadowproj, x=xshift, y=yshift)
    }
    if(info$sunaz > 0 & info$sunaz < 90){
      angle = 90-info$sunaz
      angle = 10
      yshift = round((sin((pi/180)*angle) * shiftlen *-1)/reso)*reso
      xshift = round((cos((pi/180)*angle) * shiftlen *-1)/reso)*reso
      shadowproj = shift(shadowproj, x=xshift, y=yshift)
    }
    return(shadowproj)
  }
  
  for(m in 1:length(shiftlen)){
    if(m == 1){mergeit = "r1"} else {mergeit = paste(mergeit,",r",m, sep="")}
    dothis = paste("r",m,"=shiftit(shadowproj,shiftlen[m],info,reso)", sep="")
    eval(parse(text=dothis))
    if(m == length(shiftlen)){mergeit = paste("shadowproj = mosaic(",mergeit,",fun=max,na.rm=T)", sep="")}
  }
  
  shadowproj1 = shiftit(shadowproj,shiftlen[1],info,reso)
  shadowproj2 = shiftit(shadowproj,shiftlen[2],info,reso)
  shadowproj3 = shiftit(shadowproj,shiftlen[3],info,reso)
  shadowproj4 = shiftit(shadowproj,shiftlen[4],info,reso)
  shadowproj5 = shiftit(shadowproj,shiftlen[5],info,reso)
  shadowproj6 = shiftit(shadowproj,shiftlen[6],info,reso)

  shadowproj = mosaic(shadowproj1,shadowproj2,shadowproj3,shadowproj4,shadowproj5,shadowproj6,fun=max) 
  shadowproj = extend(shadowproj, cloud, value=0)
  shadowproj = crop(shadowproj, cloud)
  
  shadow = setValues(ref,shadow)
  shadow = shadow*shadowproj
  shadow = as.matrix(shadow)
  
  #filter the aggregated cloud and shadow 
  clumps = .Call("ccl", shadow, PACKAGE = "SDMTools")  
  clumps = setValues(ref, clumps)
  fre = freq(clumps)
  these = which(fre[,2] < 10) #10
  values = fre[these,1]
  m = match(as.matrix(clumps), values) 
  these = which(is.na(m) == F)
  shadow[these] = 0
  shadow = setValues(ref, shadow)
  
  shadow = focal(shadow, w=matrix(1,5,5), fun=max)
  cloud = cloud*2
  cloudshadow = mosaic(cloud,shadow,fun=max)

  cloudshadow[is.na(ref)] = NA
  projection(cloudshadow) = set_projection(file)
  cloudshadow = as(cloudshadow, "SpatialGridDataFrame")
  if(test == F){outfile = sub("reflectance", "cloudmask", file)} else {outfile = sub("reflectance", "cloudmask", file)}
  writeGDAL(cloudshadow, outfile, drivername = "GTiff", type = "Byte", mvFlag=255, options="INTERLEAVE=BAND")
  
}
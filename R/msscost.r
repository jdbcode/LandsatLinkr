#' Convert DN values to surface reflectance
#'
#' Convert DN values to surface reflectance using the COST model with dark object subtraction
#' @param file The full path name of the *archv file
#' @import raster
#' @export


msscost = function(file, overwrite=F){
  
  #link to the equations to convert DN to TOA and TOA to SR
  #http://landsathandbook.gsfc.nasa.gov/data_prod/prog_sect11_3.html
  
  check = file_check(file,"dos_sr.tif",overwrite)
  if(check == 0){return(0)}
  
  refl = function(file, band, gain, bias, sunzenith, d, esun, dov){
    orig = raster(file, band)
    img = as.matrix(orig)   
    img = ((gain*img)+bias)-((gain*dov)+bias)
    img[img < 0] = 0
    img = (pi * img * (d^2))/(esun * cos(sunzenith))
    img = round(img * 10000)
    img = setValues(orig,img)
    return(img)
  }
  
  
  #find the dark object values  
  brightthresh = 120
  
  b1 = as.matrix(raster(file, 1))
  goodpix = which(is.na(b1) == F & b1 <= brightthresh)
  samp = matrix(nrow=100000,ncol=30)
  for(k in 1:length(samp[1,])){samp[,k] = sample(length(goodpix), 100000)}
  
  #setup the png darkobject file
  pngout = paste(substr(file,1,nchar(file)-4),"_drkobjv.png", sep="")
  if(file.exists(pngout) == T){unlink(pngout)}
  png(pngout,width=700, height=700)
  par(mfrow=c(2,2))
  
  
  #iterate through the bands finding the dark objects
  thresh = c(5,5,4,3)
  drkobjv = array(dim=4)
  for(b in 1:4){
    dn = as.matrix(raster(file, band=b))
    dnsamp = dn[goodpix]
    drkobj = array(dim=length(samp[1,]))
    for(g in 1:length(samp[1,])){
      r = dnsamp[samp[,g]]
      count = table(r)   #250
      count1 = c(count[2:length(count)],0)
      shift = (count1-count)
      valu = as.numeric(rownames(shift))
      goods = which(shift >= thresh[b])
      drkobj[g] = valu[goods[1]]
    }
    
    finaldrkobj = round(mean(drkobj))
    hist(dn, breaks=256, ylim=c(0,20000),
         main=paste("Band",b,"dark object value =", finaldrkobj), col="black", xlab="DN")
    abline(v = finaldrkobj, col = "red")
    drkobjv[b] = finaldrkobj
  }    
  dev.off()  
  
  #read in the image metadata  
  info = get_metadata(file)
 
  #define esun values for mss (chander et al 2009 summary of current radiometric calibration coefficients... RSE 113)
  if(info$sensor == "LANDSAT_1"){esun = c(1823,1559,1276,880.1)}
  if(info$sensor == "LANDSAT_2"){esun = c(1829,1539,1268,886.6)}
  if(info$sensor == "LANDSAT_3"){esun = c(1839,1555,1291,887.9)}
  if(info$sensor == "LANDSAT_4"){esun = c(1827,1569,1260,866.4)}
  if(info$sensor == "LANDSAT_5"){esun = c(1824,1570,1249,853.4)}
  
  #define the earth sun distance 
  d = eudist(info$doy)
  sunzenith = info$sunzen*(pi/180) 
  
  #apply the conversion to reflectance
  b1 = refl(file,1,info$b1gain, info$b1bias, sunzenith, d, esun[1], drkobjv[1]) 
  b2 = refl(file,2,info$b2gain, info$b2bias, sunzenith, d, esun[2], drkobjv[2])
  b3 = refl(file,3,info$b3gain, info$b3bias, sunzenith, d, esun[3], drkobjv[3])
  b4 = refl(file,4,info$b4gain, info$b4bias, sunzenith, d, esun[4], drkobjv[4])
  
  #stack the bands
  img = stack(b1,b2,b3,b4)
  
  #write out
  projection(img) = set_projection(file)
  img = as(img, "SpatialGridDataFrame") 
  outfile = sub("archv", "dos_sr", file)
  writeGDAL(img, outfile, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  
  return(1)
}
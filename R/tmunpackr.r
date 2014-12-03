#' Decompress, stack, and reproject LPSG MSS images
#'
#' Decompresses, stacks, and optionally reprojects LPGS MSS images recieved from USGS EROS as .tar.gz files
#' @param file character. full path name of the surface reflectance file
#' @param outtype coded integer designating what units are desired in the output image: 1=DN, 2=radiance, 3=surface reflectance
#' @import raster
#' @import gdalUtils
#' @export

# http://earthexplorer.usgs.gov/ Landsat CDR TM and ETM+ images
tmunpackr = function(file, proj="default", reso=60){
  #tcset = "all", "tc", "tca"
  
  #set new directories
  randomstring = paste(sample(c(0:9, letters, LETTERS), 6, replace=TRUE),collapse="")
  tempdir = file.path(dirname(file),randomstring) #temp
  year = substr(basename(file),10,13)
  
  pieces = unlist(strsplit(dirname(file), "/")) #break up the directory and unlist so the pieces can be called by index
  len = length(pieces)-1 #get the ending index for "scene"
  newpieces = paste(pieces[1:len], collapse = "/") #subset the directory pieces so the last piece is the scene
  outdir = file.path(newpieces, "images", year)
  dir.create(tempdir, recursive=T, showWarnings=F)
  dir.create(outdir, recursive=T, showWarnings=F)
  
  #decompress the image and get/set files names
  untar(file, exdir=tempdir) #decompress the file
  files = list.files(tempdir, full.names=T)
  bands = grep("band", files, value=T)
  shadow = grep("cloud_shadow_qa", files, value=T) #0 okay, 255 bad
  cloud = grep("sr_cloud_qa.tif", files, value=T) #0 okay, 255 bad
  snow = grep("sr_cloud_qa.tif", files, value=T) #0 okay, 255 bad
  fmask = grep("cfmask", files, value=T) # <= 1 okay background 255
  outbase = substr(basename(file),1,16) 
  tempstack = file.path(tempdir,paste(outbase,"_tempstack.tif",sep=""))
  tempvrt = sub("tempstack.tif", "tempmask.vrt", tempstack)
  tempmask = sub("tempstack", "tempmask", tempstack)
  projstack = sub("tempstack", "projstack", tempstack)
  projmask = sub("tempstack", "projmask", tempstack)
  finalstack = file.path(outdir,paste(outbase,"_ledaps.tif", sep=""))
  finalmask = file.path(outdir,paste(outbase,"_cloudmask.tif", sep=""))
  tcfile = file.path(outdir,paste(outbase,"_tc.tif", sep=""))
  tcafile = file.path(outdir,paste(outbase,"_tca.tif", sep=""))
  outprojfile = file.path(outdir,paste(outbase,"_proj.txt", sep=""))
  
  ref = raster(bands[1]) #set a reference raster for setting values and getting projection
  origproj = projection(ref)
  
  #stack the image bands and write out
  gdalbuildvrt(gdalfile=bands, output.vrt = tempvrt, separate=T) #, tr=c(reso,reso)
  gdal_translate(src_dataset=tempvrt, dst_dataset=tempstack, of = "GTiff", co="INTERLEAVE=BAND")
  
   #s = stack(bands)
   #origproj = projection(s)
   #s = as(s, "SpatialGridDataFrame")       
   #writeGDAL(s, tempstack, drivername = "ENVI", type = "Int16", mvFlag = -9999)
  
  #make a composite cloudmask
  s = !is.na(as.matrix(raster(shadow)))
  c = !is.na(as.matrix(raster(cloud)))
  sn = !is.na(as.matrix(raster(snow)))
  f = as.matrix(raster(fmask)) <= 1
  mask = s*c*f*sn
  mask = setValues(ref,mask)
  mask = as(mask, "SpatialGridDataFrame")        #convert the raster to SGHF so it can be written using GDAL (faster than writing it with the raster package)
  writeGDAL(mask, tempmask, drivername = "GTiff", type = "Byte", mvFlag = 255, options="INTERLEAVE=BAND")
  
  s=c=sn=f=mask=0 #clear the memory
  
  #reproject the image #need to add in writing proj file for default
  if(proj == "default"){proj = origproj}
  if(proj == "albers"){proj = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"}
    write(proj, outprojfile)
    gdalwarp(srcfile=tempstack, dstfile=projstack, 
               s_srs=origproj, t_srs=proj, of="GTiff", 
               r="bilinear", srcnodata=-9999, dstnodata=-32768, multi=T, #"near"
               tr=c(reso,reso), co="INTERLEAVE=BAND")
    
    #project the mask
    gdalwarp(srcfile=tempmask, dstfile=projmask, 
             s_srs=origproj, t_srs=proj, of="GTiff", 
             r="mode", srcnodata=255, dstnodata=255, multi=T,
             tr=c(reso,reso), co="INTERLEAVE=BAND")

  

  #trim the na rows and cols
  if(proj != "default"){
    infile = projstack
    inmask = projmask
  } else {
    infile = tempstack
    inmask = tempmask
  } 
  trim_na_rowcol(infile, finalstack, inmask, finalmask)
  
  #tasseled cap
  ref = raster(finalstack, 1)
  b1 = as.matrix(raster(finalstack, 1))
  b2 = as.matrix(raster(finalstack, 2))
  b3 = as.matrix(raster(finalstack, 3))
  b4 = as.matrix(raster(finalstack, 4))
  b5 = as.matrix(raster(finalstack, 5))
  b6 = as.matrix(raster(finalstack, 6))
  
  bcoef = c(0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303)
  gcoef = c(-0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446)
  wcoef = c(0.0315, 0.2021, 0.3102, 0.1594,-0.6806, -0.6109)
  
  tcset="all" #hardwire

  bright = (b1*bcoef[1])+(b2*bcoef[2])+(b3*bcoef[3])+(b4*bcoef[4])+(b5*bcoef[5])+(b6*bcoef[6])
  green = (b1*gcoef[1])+(b2*gcoef[2])+(b3*gcoef[3])+(b4*gcoef[4])+(b5*gcoef[5])+(b6*gcoef[6])
  if(tcset == "all" | tcset == "tc"){wet = (b1*wcoef[1])+(b2*wcoef[2])+(b3*wcoef[3])+(b4*wcoef[4])+(b5*wcoef[5])+(b6*wcoef[6])}
  
  b1=b2=b3=b4=b5=b6=0
  
  #calc tc and convert to a raster
  if(tcset == "all" | tcset == "tc"){
    tcb = setValues(ref,bright)
    tcg = setValues(ref,green)
    tcw = setValues(ref,wet)
    tc = stack(tcb,tcg,tcw)
    projection(tc) = set_projection(tcfile)
    tc = as(tc, "SpatialGridDataFrame")
    writeGDAL(tc, tcfile, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  }
  
  if(tcset == "all" | tcset == "tca"){
    #calc tc angle and convert to a raster
    tca = atan(green/bright) * (180/pi) * 100
    tca = setValues(ref,tca)
    projection(tca) = set_projection(tcafile)
    tca = as(tca, "SpatialGridDataFrame")
    writeGDAL(tca, tcafile, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  }
  
  #delete temporary files
  unlink(tempdir, recursive=T, force=T)
}

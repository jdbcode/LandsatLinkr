#' Decompress, stack, and reproject OLI SR images
#'
#' Decompress, stack, and reproject OLI SR images
#' @param file character. full path name of the surface reflectance file
#' @param proj character. PROJ.4 projection definition.
#' @param overwrite logical. True will overwrite the file if it already exists, False will skip processing if output file exists. 
#' @import raster
#' @import gdalUtils
#' @import rgdal
#' @export


oliunpackr = function(file, proj="default", overwrite=F){

  # http://earthexplorer.usgs.gov/ Landsat CDR OLI images
  
  check = file_check(file,"l8sr.tif",overwrite)
  print(check)
  if(check == 0){return(0)}
  
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
  untar(file, exdir=tempdir, tar="internal") #decompress the file
  files = list.files(tempdir, full.names=T)
  bands = grep("band", files, value=T)
  fmask = grep("cfmask.tif", files, value=T) # <= 1 okay background 255
  outbase = substr(basename(file),1,16) 
  tempstack = file.path(tempdir,paste(outbase,"_tempstack.tif",sep=""))
  tempvrt = sub("tempstack.tif", "tempstack.vrt", tempstack)
  tempmask = sub("tempstack", "tempmask", tempstack)
  projstack = sub("tempstack", "projstack", tempstack)
  projmask = sub("tempstack", "projmask", tempstack)
  finalstack = file.path(outdir,paste(outbase,"_l8sr.tif", sep=""))
  finalmask = file.path(outdir,paste(outbase,"_cloudmask.tif", sep=""))
  outprojfile = file.path(outdir,paste(outbase,"_proj.txt", sep=""))
  
  ref = raster(bands[1]) #set a reference raster for setting values and getting projection
  origproj = projection(ref)
  
  #stack the image bands and write out
  gdalbuildvrt(gdalfile=bands, output.vrt = tempvrt, separate=T) #, tr=c(reso,reso)
  gdal_translate(src_dataset=tempvrt, dst_dataset=tempstack, of = "GTiff", co="INTERLEAVE=BAND")
  
  #create the mask
  mask = as.matrix(raster(fmask))
  mask = mask <= 1
  mask = mask*1 #convert from logial to numeric
  mask = setValues(ref,mask)
  mask = as(mask, "SpatialGridDataFrame")        #convert the raster to SGHF so it can be written using GDAL (faster than writing it with the raster package)
  writeGDAL(mask, tempmask, drivername = "GTiff", type = "Byte", mvFlag = 255)
  
  
  #reproject the image #need to add in writing proj file for default
  if(proj == "default"){proj = origproj}
  write(proj, outprojfile)
  gdalwarp(srcfile=tempstack, dstfile=projstack, 
           s_srs=origproj, t_srs=proj, of="GTiff", 
           r="bilinear", srcnodata=-9999, dstnodata=-32768, multi=T, #"near"
           tr=c(30,30), co="INTERLEAVE=BAND")
  
  #project the mask
  gdalwarp(srcfile=tempmask, dstfile=projmask, 
           s_srs=origproj, t_srs=proj, of="GTiff", 
           r="mode", srcnodata=255, dstnodata=255, multi=T,
           tr=c(30,30), co="INTERLEAVE=BAND")
  
  
  #trim the na rows and cols
  trim_na_rowcol(projstack, finalstack, projmask, finalmask)
  
  #delete temporary files
  unlink(tempdir, recursive=T, force=T)
  return(1)
}


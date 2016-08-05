#' Decompress, stack, and reproject TM/ETM+ SR images
#'
#' Decompress, stack, and reproject TM/ETM+ SR images
#' @param file character. full path name of the surface reflectance file
#' @param proj character. PROJ.4 projection definition.
#' @param overwrite logical. True will overwrite the file if it already exists, False will skip processing if output file exists. 
#' @import raster
#' @import gdalUtils
#' @import rgdal
#' @export


tmunpackr = function(file, proj="default", overwrite=F){
  # http://earthexplorer.usgs.gov/ Landsat CDR TM and ETM+ images
  
  check = file_check(file,"ledaps.tif",overwrite)
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
  shadow = grep("cloud_shadow_qa.tif", files, value=T) #0 okay, 255 bad
  cloud = grep("sr_cloud_qa.tif", files, value=T) #0 okay, 255 bad
  snow = grep("sr_snow_qa.tif", files, value=T) #0 okay, 255 bad
  fmask = grep("cfmask.tif", files, value=T) # <= 1 okay background 255
  outbase = substr(basename(file),1,16) 
  tempstack = file.path(tempdir,paste(outbase,"_tempstack.tif",sep=""))
  tempvrt = sub("tempstack.tif", "tempstack.vrt", tempstack)
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
  
  #make a composite cloudmask
  s = as.matrix(raster(shadow))
  c = as.matrix(raster(cloud))
  sn = as.matrix(raster(snow))
  f = as.matrix(raster(fmask))
  
  check = s[1,1] # if is.na(check) == T new else old
  if(is.na(check) == T){s = is.na(s)} else {s = !is.na(s)}
  check = c[1,1] # if is.na(check) == T new else old
  if(is.na(check) == T){c = is.na(c)} else {c = !is.na(c)}
  check = sn[1,1] # if is.na(check) == T new else old
  if(is.na(check) == T){sn = is.na(sn)} else {sn = !is.na(sn)}
  f = f <= 1
  mask = s*c*f*sn
  mask = setValues(ref,mask)
  mask = as(mask, "SpatialGridDataFrame")        #convert the raster to SGHF so it can be written using GDAL (faster than writing it with the raster package)
  writeGDAL(mask, tempmask, drivername = "GTiff", type = "Byte", mvFlag = 255, options="INTERLEAVE=BAND")
  
  s=c=sn=f=mask=0 #clear the memory
  
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
  
  
  bright = (b1*bcoef[1])+(b2*bcoef[2])+(b3*bcoef[3])+(b4*bcoef[4])+(b5*bcoef[5])+(b6*bcoef[6])
  green = (b1*gcoef[1])+(b2*gcoef[2])+(b3*gcoef[3])+(b4*gcoef[4])+(b5*gcoef[5])+(b6*gcoef[6])
  wet = (b1*wcoef[1])+(b2*wcoef[2])+(b3*wcoef[3])+(b4*wcoef[4])+(b5*wcoef[5])+(b6*wcoef[6])

  b1=b2=b3=b4=b5=b6=0

  tcb = matrix_to_raster(finalstack,bright)
  tcg = matrix_to_raster(finalstack,green)
  tcw = matrix_to_raster(finalstack,wet)
  wet=0  
  tca = atan(green/bright) * (180/pi) * 100
  bright=green=0  
  tca = matrix_to_raster(finalstack,tca)
  
  temptcb = file.path(tempdir,paste(outbase,"_temptcb.tif",sep=""))
  temptcg = file.path(tempdir,paste(outbase,"_temptcg.tif",sep=""))
  temptcw = file.path(tempdir,paste(outbase,"_temptcw.tif",sep=""))
  projection(tcb) = set_projection(tcfile)
  projection(tcg) = set_projection(tcfile)
  projection(tcw) = set_projection(tcfile)
  tc = as(tcb, "SpatialGridDataFrame")
  tcb=0
  writeGDAL(tc, temptcb, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  tc = as(tcg, "SpatialGridDataFrame")
  tcg=0
  writeGDAL(tc, temptcg, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  tc = as(tcw, "SpatialGridDataFrame")
  tcw=0
  writeGDAL(tc, temptcw, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  
  bands = c(temptcb,temptcg,temptcw)
  temptcs = file.path(tempdir,paste(outbase,"_temptcstack.vrt",sep=""))
  gdalbuildvrt(gdalfile=bands, output.vrt = temptcs, separate=T) #, tr=c(reso,reso)
  gdal_translate(src_dataset=temptcs, dst_dataset=tcfile, of = "GTiff", co="INTERLEAVE=BAND")
  
  projection(tca) = set_projection(tcfile)
  tc = as(tca, "SpatialGridDataFrame")
  tca=0
  writeGDAL(tc, tcafile, drivername = "GTiff", type = "Int16", mvFlag = -32768, options="INTERLEAVE=BAND")
  
  #delete temporary files
  unlink(tempdir, recursive=T, force=T)
}


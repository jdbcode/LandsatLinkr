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

  #file = "D:/work/proj/llr_dev/collection1/oli/wrs2/044034/targz/LC080440342018020901T1-SC20180309132638.tar.gz"
  #proj= "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  
  
  # http://earthexplorer.usgs.gov/ Landsat CDR OLI images
  
  check = file_check(file,"l8sr.tif",overwrite)
  if(check == 0){return(0)}
  
  #make a random work folder
  randomstring = paste(sample(c(0:9, letters, LETTERS), 6, replace=TRUE),collapse="")
  tempdir = file.path(dirname(file),randomstring) #temp
  dir.create(tempdir, recursive=T, showWarnings=F)
  
  #decompress the image 
  untar(file, exdir=tempdir, tar="internal") #decompress the file
  files = list.files(tempdir, full.names=T)
  
  #are there files in the decompressed archive
  if(length(files) == 0){
    unlink(tempdir, recursive=T, force=T)
    stop(paste('There were no files found in the decompressed archive:', file))
  }
  
  #are these sr files
  filesBname = basename(files)
  srFilesTest = grep("_sr_", filesBname, value=T)
  if(length(srFilesTest) == 0){
    unlink(tempdir, recursive=T, force=T)
    stop(paste('There were no USGS ESPA Surface Reflectance files found in the decompressed archive:', file))
  }
  
  #is this collection 1 or pre-collection
  name = filesBname[1]
  isCollection = substr(filesBname[1], 5,5) == '_'
  
  #get the sr bands
  bands = sort(grep("band", files, value=T))
  
  if(isCollection){
    #new version: LXSS_LLLL_PPPRRR_YYYYMMDD_yyyymmdd_CX_TX_prod_band.ext
    year = substr(name,18,21)
    pieces = unlist(strsplit(dirname(file), "/")) #break up the directory and unlist so the pieces can be called by index
    len = length(pieces)-1 #get the ending index for "scene"
    newpieces = paste(pieces[1:len], collapse = "/") #subset the directory pieces so the last piece is the scene
    outdir = file.path(newpieces, "images", year)
    dir.create(outdir, recursive=T, showWarnings=F)
    
    # create the mask
    pixelqafile = grep("pixel_qa.tif", files, value=T)
    pixelqar = getValues(raster(pixelqafile))
    #mask = as.numeric(pixelqar == 322 | pixelqar == 386 | pixelqar == 324 | pixelqar == 388 | pixelqar == 836 | pixelqar == 900)
    shadow = bitwAnd(pixelqar, 8) == 0  # shadow
    snow = bitwAnd(pixelqar, 16) == 0  # snow
    clouds = bitwAnd(pixelqar, 32) == 0 # clouds
    mask = shadow*snow*clouds

    clouds = snow = shadow = pixelqar = 0 #clear memory
    
    # make the basename for final output files
    mtlfile = grep("MTL.txt", files, value=T)
    tbl = unlist(read.delim(mtlfile, header=F, skipNul=T))
    string = as.character(grep("LANDSAT_SCENE_ID = ", tbl, value=T))
    pieces = unlist(strsplit(string, " "))
    sceneid = pieces[length(pieces)]
    outbase = substr(sceneid,1,16)
    
  } else{
    #name = 'LC80380292015214LGN00.xml'
    outbase = substr(name,1,16) # need to get this from the filename, because pre-collection does not include an mtl file
    year = substr(name,10,13)
    
    pieces = unlist(strsplit(dirname(file), "/")) #break up the directory and unlist so the pieces can be called by index
    len = length(pieces)-1 #get the ending index for "scene"
    newpieces = paste(pieces[1:len], collapse = "/") #subset the directory pieces so the last piece is the scene
    outdir = file.path(newpieces, "images", year)
    dir.create(outdir, recursive=T, showWarnings=F)
    
    # create the mask
    fmask = grep("cfmask.tif", files, value=T) # <= 1 okay background 255
    mask = as.matrix(raster(fmask))
    mask = mask <= 1
    mask = mask*1 #convert from logial to numeric
  }
  
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
  
  mask = setValues(ref,mask)
  mask = as(mask, "SpatialGridDataFrame")        #convert the raster to SGHF so it can be written using GDAL (faster than writing it with the raster package)
  writeGDAL(mask, tempmask, drivername = "GTiff", type = "Byte", mvFlag = 255)
  mask = 0 # clear memory
  
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


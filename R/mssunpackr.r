#' Decompress, stack, and reproject LPSG MSS images
#'
#' Decompresses, stacks, and optionally reprojects LPGS MSS images recieved from USGS EROS as .tar.gz files
#' @param file Filename of LPGS Landsat MSS image filename (full system path to file) 
#' @param proj PROJ.4 projection definition. By default no projection will take place. Optionally specify a CRS projection string or "albers" for the USGS version of Albers Equal Area Conic 
#' @param reso numeric. the target pixel size for the output image
#' @import raster
#' @import gdalUtils
#' @export


mssunpackr = function(file, proj, reso=60){
  
  randomstring = paste(sample(c(0:9, letters, LETTERS), 6, replace=TRUE),collapse="")
  tempdir = file.path(dirname(file),randomstring) #temp
  untar(file, exdir=tempdir) #decompress the file
  mtlfile = list.files(tempdir, pattern = "MTL.txt", full.names = T, recursive = T)
  tbl = unlist(read.delim(mtlfile, header=F, skipNul=T))
  dtype = as.character(grep("DATA_TYPE = ", tbl, value=T))
  l1ttest = length(grep("L1T", dtype))
  if(l1ttest == 1){
    #find files and make new file names and directories
    allfiles = list.files(tempdir, full.names=T) #find all the decompressed files 
    tiffiles = allfiles[grep("TIF",allfiles)] #subset the tif image files
    verfile = allfiles[grep("VER.txt",allfiles)]
    otherfiles = allfiles[grep("TIF",allfiles, invert=T)] #subset the other files
    filebase = basename(tiffiles[1]) #get the basename
    filedir = dirname(file) #get the directory
    year = substr(filebase, 10, 13) #get the year
    pieces = unlist(strsplit(filedir, "/")) #break up the directory and unlist so the pieces can be called by index
    len = length(pieces)-1 #get the ending index for "scene"
    newpieces = paste(pieces[1:len], collapse = "/") #subset the directory pieces so the last piece is the scene
    imgid = substr(filebase, 1, 16)
    name = paste(imgid, "_archv.tif", sep = "") #define the new file basename 
    newdir = file.path(newpieces, "images", year, fsep = .Platform$file.sep) #define the new directory
    
    #make output filenames
    tempstack = file.path(tempdir,paste(imgid, "_tempstack.tif", sep = ""))
    projstack = sub("tempstack", "projstack", tempstack)
    finalstack = file.path(newdir, name) #define the new full filename of the output image 
    outprojfile = sub("archv.tif", "proj.txt", finalstack)
    
    #deal with the ancillary file names
    junk = substr(filebase, 17,21)
    baseotherfiles = basename(otherfiles)
    for(h in 1:length(otherfiles)){baseotherfiles[h] = sub(junk, "", baseotherfiles[h])}
    newotherfiles =  file.path(newdir, baseotherfiles, fsep = .Platform$file.sep) #define the new filenames for associated files
    
    ref = raster(tiffiles[1]) #load a file to get original projection
    s = stack(tiffiles[1],tiffiles[2],tiffiles[3],tiffiles[4])
    img = as.array(s)    
    b1bads = img[,,1]>1 #!=0
    b2bads = img[,,2]>1 #!=0
    b3bads = img[,,3]>1 #!=0
    b4bads = img[,,4]>1 #!=0
    bads = b1bads*b2bads*b3bads*b4bads
    
    if(length(which(bads==0))/ncell(ref) > 0.75){
      delete_files(file, 5)
      unlink(c(tempdir), recursive=T, force=T)
      return(0)
    }
    
    img[,,1] = img[,,1]*bads
    img[,,2] = img[,,2]*bads
    img[,,3] = img[,,3]*bads
    img[,,4] = img[,,4]*bads

    dir.create(newdir, recursive=T, showWarnings=F) #make a new output directory
    cloudtest = img[,,1] > 130
    percent = round((length(which(cloudtest == T))/length(which((img[,,1]>0) == T)))*100, digits=4)
    tbl = unlist(read.delim(verfile, header=F, skipNul=T))
    rmseline = as.character(grep("Scene RMSE: ", tbl, value=T))
    rmse = as.numeric(unlist(strsplit(rmseline, " "))[3])
    outfile = sub("archv.tif", "cloud_rmse.csv", finalstack)
    write(c(finalstack,percent,rmse),file=outfile, ncolumns=3, sep=",")
    
    s = setValues(s,img)    
    origproj = projection(s)
    s = as(s, "SpatialGridDataFrame")       
    writeGDAL(s, tempstack, drivername = "GTiff", options="INTERLEAVE=BAND", type = "Byte", mvFlag = 0) #, drivername = "GTiff"
    
    #reproject the image #need to add in writing proj file for default
    #if(proj != "default"){
      #if(proj == "albers"){proj = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"}
      write(proj, outprojfile)
      gdalwarp(srcfile=tempstack, dstfile=projstack, 
               s_srs=origproj, t_srs=proj, of="GTiff",
               r="near", srcnodata=0, dstnodata=0, multi=T,
               tr=c(reso,reso), co="INTERLEAVE=BAND")
#     } else {
#       write(origproj, outprojfile)
#       gdalwarp(srcfile=tempstack, dstfile=projstack, 
#                s_srs=origproj, t_srs=origproj, of="GTiff",
#                r="near", srcnodata=0, dstnodata=0, multi=T,
#                tr=c(reso,reso), co="INTERLEAVE=BAND")
#     }
    
    #trim the na rows and cols
    if(proj != "default"){infile = projstack} else {infile = tempstack} 
    trim_na_rowcol(infile, finalstack, "null", "null") 
    
    file.rename(otherfiles,newotherfiles) #move the associated files
    unlink(tempdir, recursive=T, force=T) #delete the temp directory
    return(1)
  } else {
    #outfile = file.path(dirname(file), paste(imgid,"_targz_delete.csv",sep=""))
    #write(c(imgid,"L1G"),file=outfile, ncolumns=2, sep=",")
    delete_files(file,1)
    unlink(c(tempdir), recursive=T, force=T)
    return(0)
  }
}
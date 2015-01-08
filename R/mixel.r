#' Composite images 
#'
#' Composite images
#' @param msswrs1dir character. list of mss wrs1 directory paths
#' @param msswrs2dir character. list of mss wrs2 directory paths
#' @param tmwrs2dir character. list of tm wrs2 directory path
#' @param index character. spectral index to make composites for. options: "tca", "tcb", "tcg", "tcw"
#' @param outdir character. path to output directory
#' @param runname character. unique name for the composite set
#' @param useareafile character. path to usearea file
#' @param doyears ??? what years to composite
#' @param order character. how to order the images options "sensor_and_doy" and "doy"
#' @param overlap character. how to deal with overlapping images. options: "mean"
#' @import raster
#' @import gdalUtils
#' @export


mixel = function(msswrs1dir,msswrs2dir,tmwrs2dir,index,outdir,runname,useareafile,doyears="all",order="sensor_and_doy",overlap="mean"){
  #year = "all" or any year that you want
  
  mixel_find = function(files, refimg){
    
    #get the extents of the files
    info = matrix(ncol = 4, nrow=length(files))
    print("Getting image extents")
    for(i in 1:length(files)){ 
      print(i)
      img = raster(files[i])
      ext = extent(img)
      info[i,1] = ext@xmin
      info[i,2] = ext@xmax
      info[i,3] = ext@ymin
      info[i,4] = ext@ymax
    }
    
    text = extent(raster(refimg))  
    these = which(info[,3] < text@ymax & info[,4] > text@ymin & info[,2] > text@xmin & info[,1] < text@xmax) 
    goods = files[these]
    return(goods)
  }
  
  
  mixel_mask = function(imgfile, search, refimg){
    maskfile = sub(search, "cloudmask.tif",imgfile)
    print(maskfile)
    print(file.exists(maskfile))
    img = raster(imgfile)
    mask = raster(maskfile)
    
    imgex = alignExtent(img, refimg, snap="near")
    maskex = alignExtent(mask, refimg, snap="near")
    extent(img) = imgex
    extent(mask) = maskex
    
    overlap = intersect(img, mask)
    mask = crop(mask, overlap)
    img = crop(img, overlap)
    
    bads = which(values(mask) == 0)
    img[bads] = NA
    return(img)
  }
  
  
  if(index == "tca"){search="tca.tif"}
  if(index == "tcb"){search="tc.tif"}
  if(index == "tcg"){search="tc.tif"}
  if(index == "tcw"){search="tc.tif"}
  
  msswrs1imgdir = file.path(msswrs1dir[i],"images")
  msswrs2imgdir = file.path(msswrs2dir,"images")
  tmwrs2imgdir = file.path(tmwrs2dir,"images")
  
  for(i in 1:length(msswrs1dir)){
    if(i == 1){msswrs1files = list.files(msswrs1imgdir[i], search, recursive=T, full.names=T)} else {
      msswrs1files = c(msswrs1files,list.files(msswrs1imgdir[i], search, recursive=T, full.names=T))
    }
  }
  for(i in 1:length(msswrs2dir)){
    if(i == 1){msswrs2files = list.files(msswrs2imgdir[i], search, recursive=T, full.names=T)} else {
      msswrs2files = c(msswrs2files,list.files(msswrs2imgdir[i], search, recursive=T, full.names=T))
    }
  }  
  for(i in 1:length(msswrs2dir)){
    if(i == 1){tmwrs2files = list.files(tmwrs2imgdir[i], search, recursive=T, full.names=T)} else {
      tmwrs2files = c(tmwrs2files, list.files(tmwrs2imgdir[i], search, recursive=T, full.names=T))
    }
  }
  
  files = c(msswrs1files,msswrs2files,tmwrs2files)
  files = mixel_find(files, useareafile)
  
  #create an output directory
  dir.create(outdir, recursive=T, showWarnings=F)
  
  #extract from info from the filenames
  filebase = basename(files)
  years = substr(filebase, 10, 13)
  days = substr(filebase, 14,16)
  sensor = substr(filebase, 1,3)
  yearsort = sort(unique(years))
  medday = median(as.numeric(days))
  
  #load in the reference image (tsa usearea file)
  refimg = raster(useareafile)
  refimg[values(refimg) == 1] = NA
  
  if(doyears == "all"){uni=yearsort} else {
    theseyears = match(doyears,yearsort)
    uni = yearsort[theseyears]
  }
  
  #for all the unique year make a composite
  for(i in 1:length(uni)){
    ptm <- proc.time()
    if(is.na(uni[i] == T)){next}
    print(paste("working on year:", uni[i]))
    these = which(years == uni[i])
    theseimgs = files[these]
    thesedays = days[these]
    thesesensors = sensor[these]
    meddif = abs(as.numeric(thesedays)-medday)
    
    #day of year order - no sensor consideration
    if(order == "doy"){
      difsort = sort(meddif, index.return = T)
      imgorder = theseimgs[difsort$ix]
    }
    
    #use this section if TM is to always come first in the img merge order
    if(order == "sensor_and_doy"){
      #separate the sensors
      tmid = which(thesesensors == "LT5" | thesesensors == "LT4")
      etmid = which(thesesensors == "LE7")
      mssid = which(thesesensors != "LE7" & thesesensors != "LT5" & thesesensors != "LT4")
      
      #extract images for each sensor
      tmimg = theseimgs[tmid]
      etmimg = theseimgs[etmid]
      mssimg = theseimgs[mssid]
      
      #sort each sensor by image date
      tmmeddiff = sort(meddif[tmid], index.return = T)
      etmmeddiff = sort(meddif[etmid], index.return = T)
      mssmeddiff = sort(meddif[mssid], index.return = T)
      
      #reorder the images
      tmorder = tmimg[tmmeddiff$ix]
      etmorder = etmimg[etmmeddiff$ix]
      mssorder = mssimg[mssmeddiff$ix]
      
      #out them all together
      imgorder = c(tmorder, etmorder, mssorder)
    }
    len = length(these)
    print("merging files:")
    for(m in 1:len){
      print(basename(imgorder[m]))
      if(m == 1){mergeit = "r1"} else {mergeit = paste(mergeit,",r",m, sep="")}
      #run the image prep function on-the-fly
      dothis = paste("r",m,"=mixel_mask(imgorder[",m,"], search, refimg)", sep="")
      eval(parse(text=dothis))
      if(m == len){
        if(overlap == "order"){mergeit = paste("newimg = merge(",mergeit,")", sep="")} #,refimg
        if(overlap == "mean"){mergeit = paste("newimg = mosaic(",mergeit,",fun=mean,na.rm=T)", sep="")}  #refimg
        if(overlap == "median"){mergeit = paste("newimg = mosaic(",mergeit,",fun=median,na.rm=T)", sep="")}
      }
    }
    
    #run the merge function
    if(len == 1){newimg = r1} else {eval(parse(text=mergeit))} #only run merge it if there are multiple files to merge
    
    #name the new file
    newbase = paste(uni[i],"_",runname,"_",index,"_composite.bsq", sep="")  #.tif #tsa,"_",
    outimgfile = file.path(outdir,newbase)  
    outtxtfile = sub("composite.bsq", "composite_img_list.csv", outimgfile)
    imgorder = data.frame(imgorder)
    colnames(imgorder) = "File"
    write.csv(imgorder, file=outtxtfile)
    
    #crop it and set na values to 0
    newimg = round(crop(newimg, refimg))
    newimg = extend(newimg, refimg, value=0)
    newimg[(values(refimg) == 0)] = 0
    newimg[is.na(newimg)] = 0

    #write out the new image
    projection(newimg) = set_projection(files[1])
    newimg = as(newimg, "SpatialGridDataFrame") #convert the raster to SGHF so it can be written using GDAL (faster than writing it with the raster package)
    writeGDAL(newimg, outimgfile, drivername = "ENVI", type = "Int16", mvFlag = -32768) #, options="INTERLEAVE=BAND"
    
    print(proc.time()-ptm) 
  }
  
  bname = paste(runname,"_",index,"_composite.bsq", sep="")
  bands = sort(list.files(outdir, bname, full.names=T))
  fullnametif = file.path(outdir,bname)
  fullnamevrt = sub(".bsq", ".vrt", fullnametif)
  gdalbuildvrt(gdalfile=bands, output.vrt = fullnamevrt, separate=T) #, tr=c(reso,reso)
  gdal_translate(src_dataset=fullnamevrt, dst_dataset=fullnametif, of = "ENVI") #, co="INTERLEAVE=BAND"
  unlink(fullnamevrt)
  
}

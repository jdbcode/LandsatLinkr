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
#' @import ggplot2
#' @import plyr
#' @export


mixel = function(msswrs1dir,msswrs2dir,tmwrs2dir,oliwrs2dir,index,outdir,runname,useareafile,doyears="all",order="none",overlap="mean"){

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
  
  mixel_mask = function(imgfile, refimg, index){ #search,
    if(index == "tca" | index == "tcb"){band=1}
    if(index == "tcg"){band=2}
    if(index == "tcw"){band=3}
    
    sensor = substr(basename(imgfile), 1,2)
    if(sensor == "LM"){maskbit = "_cloudmask_30m.tif"} else {maskbit = "_cloudmask.tif"}
    
    maskfile = file.path(dirname(imgfile), paste(substr(basename(imgfile),1,16),maskbit,sep=""))
    print(basename(maskfile))
    print(file.exists(maskfile))
    img = raster(imgfile, band=band)
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
  
  mixel_composite = function(outdir, files, runname, index, doyears, order, useareafile, overlap, offset, adj=NULL, offsetrun){
    
    #extract some info from the filenames
    filebase = basename(files)
    years = substr(filebase, 10, 13)
    days = substr(filebase, 14,16)
    sensor = substr(filebase, 1,3)
    yearsort = sort(unique(years))
    medday = median(as.numeric(days))
    
    #load in the reference image (usearea file)
    refimg = raster(useareafile)
    newimgnas = values(refimg) == 0
    #nas = values(refimg) == 1
    
    if(doyears == "all"){uni=yearsort} else {
      theseyears = match(doyears,yearsort)
      uni = yearsort[theseyears]
    }
    
    #for all the unique year make a composite
    for(i in 1:length(uni)){
      ptm <- proc.time()
      if(is.na(uni[i] == T)){next}
      print(paste("working on year:", uni[i]))
      
      refimg = raster(useareafile)
      #refimg[nas] = NA
      
      these = which(years == uni[i])
      theseimgs = files[these]
      thesedays = days[these]
      thesesensors = sensor[these]
      meddif = abs(as.numeric(thesedays)-medday)
      
      if(order == "none"){imgorder = theseimgs}
      
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
        dothis = paste("r",m,"=mixel_mask(imgorder[",m,"], refimg, index)", sep="") 
        eval(parse(text=dothis))
        if(m == len){
          if(overlap == "order"){mergeit = paste("newimg = merge(",mergeit,")", sep="")}
          if(overlap == "mean"){mergeit = paste("newimg = mosaic(",mergeit,",fun=mean,na.rm=T,tolerance=0.5)", sep="")}
          if(overlap == "median"){mergeit = paste("newimg = mosaic(",mergeit,",fun=median,na.rm=T)", sep="")}
        }
      }
      
      #run the merge function
      if(len == 1){newimg = r1} else {eval(parse(text=mergeit))} #only run merge it if there are multiple files to merge
      
      #name the new file
      newbase = paste(uni[i],"_",runname,"_",index,"_composite.bsq", sep="")
      outimgfile = file.path(outdir,newbase)  
      outtxtfile = sub("composite.bsq", "composite_img_list.csv", outimgfile)
      imgorder = data.frame(imgorder)
      colnames(imgorder) = "File"
      write.csv(imgorder, file=outtxtfile)
      
      newimg = round(crop(newimg, refimg))
      newimg = extend(newimg, refimg, value=NA)
      newimg[newimgnas] = NA
      if(is.null(adj) == F){newimg = newimg + raster(adj)}
      if(offsetrun == F){
        newimg[newimgnas] = 0
        newimg[is.na(newimg)] = 0
      }
      
      #write out the new image
      projection(newimg) = set_projection(files[1])
      writeRaster(newimg, outimgfile, format="ENVI", datatype = "INT2S",overwrite=T)
      envifilename = sub("bsq","envi",outimgfile)
      envixmlfile = paste(envifilename,".aux.xml",sep="")
      bsqxmlfile = sub("envi","bsq",envixmlfile)
      file.rename(envifilename,outimgfile)
      file.rename(envixmlfile,bsqxmlfile)
      
      #clean the temp directory      
      tempdir = dirname(rasterTmpFile())
      #print(paste("this is the tempdir:",tempdir))      
      tempfiles = list.files(tempdir,full.names=T)
      #print(paste("there are:",length(tempfiles),"files found"))
      unlink(tempfiles)
      tempfiles = list.files(tempdir,full.names=T)
      #print(paste("after attempting to delete them, there are:",length(tempfiles),"files"))
      print(proc.time()-ptm) 
    }
  }
  
  pixel_level_offset = function(sensor, dep_files, ref_files){ 
    #start offset finding procedure
    
    if(sensor == "mss"){offsetdir = file.path(outdir,"mss_offset");deprunname="lm"; refrunname="lt"; meandiffilebname = "mss_mean_dif.bsq"}
    if(sensor == "oli"){offsetdir = file.path(outdir,"oli_offset");deprunname="lc"; refrunname="le"; meandiffilebname = "oli_mean_dif.bsq"}
    dir.create(offsetdir, recursive=T, showWarnings=F)
    #composite mss calibration images
    mixel_composite(offsetdir, dep_files, runname=deprunname,index=index, doyears=doyears, order=order, useareafile=useareafile, overlap=overlap, offsetrun=T)
    #composite tm calibration images
    mixel_composite(offsetdir, ref_files, runname=refrunname,index=index, doyears=doyears, order=order, useareafile=useareafile, overlap=overlap, offsetrun=T)
    
    #find the composite images, check and sort them
    offsetfiles = list.files(offsetdir, ".bsq$", full.names=T)
    
    #separate the offset files by dependent and reference
    depoffsetfiles = offsetfiles[which(substr(basename(offsetfiles), 6,7) == deprunname)]
    refoffsetfiles = offsetfiles[which(substr(basename(offsetfiles), 6,7) == refrunname)]
    #check to make sure there is a matching number of dependent and reference files and order them
    if(length(depoffsetfiles) == length(refoffsetfiles)){
      depyear = substr(basename(depoffsetfiles),1,4)
      refyear = substr(basename(refoffsetfiles),1,4)
      depoffsetfiles = depoffsetfiles[order(depyear)]
      refoffsetfiles = refoffsetfiles[order(refyear)]
    } else {stop("there is not a matching number of dependent and reference years in the offset folder, can't continue")}
    
    #find the mean pixel-wise difference between dependent and reference images
    for(i in 1:length(depoffsetfiles)){
      if(i == 1){dif = raster(refoffsetfiles[i]) - raster(depoffsetfiles[i])} else{
        dif = sum(dif,(raster(refoffsetfiles[i]) - raster(depoffsetfiles[i])), na.rm=T)
      }
    }
    meandiforig = round(dif/length(depoffsetfiles))
    
    #write out the mean pixel-wise difference file
    projection(meandiforig) = set_projection(files[1])
    meandiffile = file.path(offsetdir,meandiffilebname)
    writeRaster(meandiforig, meandiffile, format="ENVI", datatype = "INT2S",overwrite=T)
    envifilename = sub("bsq","envi",meandiffile)
    envixmlfile = paste(envifilename,".aux.xml",sep="")
    bsqxmlfile = sub("envi","bsq",envixmlfile)
    file.rename(envifilename,meandiffile)
    file.rename(envixmlfile,bsqxmlfile)
    
    meandif=0 #memory
    
    #     #adjust the dependent composites to reflect the offset
    #     for(i in 1:length(depoffsetfiles)){
    #       print(paste("adjusting dependent file:", depoffsetfiles[i]))
    #       r = raster(depoffsetfiles[i]) + meandiforig #make the adjustment
    #       
    #       #write out the file
    #       projection(r) = set_projection(files[1])
    #       outimgfile = sub("composite.bsq","composite_adj.bsq",depoffsetfiles[i])
    #       writeRaster(r, outimgfile, format="ENVI", datatype = "INT2S",overwrite=T)
    #       envifilename = sub("bsq","envi",outimgfile)
    #       envixmlfile = paste(envifilename,".aux.xml",sep="")
    #       bsqxmlfile = sub("envi","bsq",envixmlfile)
    #       file.rename(envifilename,outimgfile)
    #       file.rename(envixmlfile,bsqxmlfile)
    #     }
    return(meandiffile)
  }
  
  print(paste("working on index:",index))
  
  #create some search terms depending on index
  if(index == "tca"){msssearch="tca_30m.tif$"; tmsearch="tca.tif$"; olisearch="tca.tif$"}
  if(index == "tcb"){msssearch="tc_30m.tif$"; tmsearch="tc.tif$"; olisearch="tc.tif$"}
  if(index == "tcg"){msssearch="tc_30m.tif$"; tmsearch="tc.tif$"; olisearch="tc.tif$"}
  if(index == "tcw"){msssearch="tc_30m.tif$"; tmsearch="tc.tif$"; olisearch="tc.tif$"}
  
  #add "images" to scene id to ensure we only get matches from the "images" directory
  msswrs1imgdir = file.path(msswrs1dir,"images")
  msswrs2imgdir = file.path(msswrs2dir,"images")
  tmwrs2imgdir = file.path(tmwrs2dir,"images")
  oliwrs2imgdir = file.path(oliwrs2dir,"images")
  
  #find all the files from all the given directories
  for(i in 1:length(msswrs1dir)){
    if(i == 1){msswrs1files = list.files(msswrs1imgdir[i], msssearch, recursive=T, full.names=T)} else {
      msswrs1files = c(msswrs1files,list.files(msswrs1imgdir[i], msssearch, recursive=T, full.names=T))
    }
  }
  for(i in 1:length(msswrs2dir)){
    if(i == 1){msswrs2files = list.files(msswrs2imgdir[i], msssearch, recursive=T, full.names=T)} else {
      msswrs2files = c(msswrs2files,list.files(msswrs2imgdir[i], msssearch, recursive=T, full.names=T))
    }
  }  
  for(i in 1:length(tmwrs2dir)){
    if(i == 1){tmwrs2files = list.files(tmwrs2imgdir[i], tmsearch, recursive=T, full.names=T)} else {
      tmwrs2files = c(tmwrs2files, list.files(tmwrs2imgdir[i], tmsearch, recursive=T, full.names=T))
    }
  }
  for(i in 1:length(oliwrs2dir)){
    if(i == 1){oliwrs2files = list.files(oliwrs2imgdir[i], tmsearch, recursive=T, full.names=T)} else {
      oliwrs2files = c(oliwrs2files, list.files(oliwrs2imgdir[i], tmsearch, recursive=T, full.names=T))
    }
  }
  
  #put all the files together in a vector and check to make sure the files intersect the useareafile, if not they will be excluded
  files = c(msswrs1files,msswrs2files,tmwrs2files,oliwrs2files)
  files = mixel_find(files, useareafile)
  if(length(files)==0){stop("There were no files in the given directories that intersect the provided 'usearea file'.
                            Make sure that you provided the correct file, checked that it actually overlaps the scenes you 
                            specified for compositing, and that the projection is the same as the images.")}
  
  
  #organize files by sensor
  sensor = substr(basename(files), 1,2)
  mssfiles = files[which(sensor == "LM")]
  tmfiles = files[which(sensor == "LT")]
  etmfiles = files[which(sensor == "LE")]
  olifiles = files[which(sensor == "LC")]
  tmetmfiles = c(tmfiles,etmfiles)
  
  #find overlapping MSS and TM
  mssyears = substr(basename(mssfiles), 10, 13)
  tmyears = substr(basename(tmfiles), 10, 13)
  thesetm = which(tmyears %in% mssyears)
  overlaptmfiles = tmfiles[thesetm]
  thesemss = which(mssyears %in% tmyears)
  overlapmssfiles = mssfiles[thesemss]
  
  #find overlapping ETM+ and OLI
  oliyears = substr(basename(olifiles), 10, 13)
  etmyears = substr(basename(etmfiles), 10, 13)
  theseetm = which(etmyears %in% oliyears)
  overlapetmfiles = etmfiles[theseetm]
  theseoli = which(oliyears %in% etmyears)
  overlapolifiles = olifiles[theseoli]
  
  #make composites
  mssdir = file.path(outdir,"mss")
  olidir = file.path(outdir,"oli")
  tmdir = file.path(outdir,"tm")
  
  if(sum(length(overlapmssfiles),length(overlapolifiles)) == 0){ #if there are no overlapping MSS to TM and ETM+ to OLI then just composite without making pixel level offset adjustment
    mixel_composite(outdir, files, runname=runname,index=index, doyears=doyears, order=order, useareafile=useareafile, overlap=overlap, offsetrun=F)
  }
  if(length(overlapmssfiles) != 0){
    #get the MSS to TM mean pixel-level offset
    offsetfile = pixel_level_offset("mss", overlapmssfiles, overlaptmfiles)
    
    #make final mss composites
    print("compositing MSS data using offset adjustment")
    dir.create(mssdir, recursive=T, showWarnings=F)
    mixel_composite(mssdir, mssfiles, runname=runname,index=index, doyears=doyears, order=order, useareafile=useareafile, overlap=overlap, adj=offsetfile, offsetrun=F)
  }
  if(length(overlapolifiles) != 0){
    #get the OLI to ETM+ mean pixel-level offset
    offsetfile = pixel_level_offset("oli", overlapolifiles, overlapetmfiles)
    
    #make final mss composites
    print("compositing OLI data using offset adjustment")
    dir.create(olidir, recursive=T, showWarnings=F)
    mixel_composite(olidir, olifiles, runname=runname,index=index, doyears=doyears, order=order, useareafile=useareafile, overlap=overlap, adj=offsetfile, offsetrun=F)
  }
  if(length(tmetmfiles) != 0){
    #make final mss composites
    print("compositing TM and ETM+ data")
    dir.create(tmdir, recursive=T, showWarnings=F)
    mixel_composite(tmdir, tmetmfiles, runname=runname,index=index, doyears=doyears, order=order, useareafile=useareafile, overlap=overlap, offsetrun=F)
  }
  if(length(mssfiles) != 0 & length(overlapmssfiles) == 0){
    #make final mss composites
    print("compositing MSS data")
    dir.create(mssdir, recursive=T, showWarnings=F)
    mixel_composite(mssdir, mssfiles, runname=runname,index=index, doyears=doyears, order=order, useareafile=useareafile, overlap=overlap, offsetrun=F)
  }
  if(length(olifiles) != 0 & length(overlapolifiles) == 0){
    #make final mss composites
    print("compositing OLI data")
    dir.create(olidir, recursive=T, showWarnings=F)
    mixel_composite(olidir, olifiles, runname=runname,index=index, doyears=doyears, order=order, useareafile=useareafile, overlap=overlap, offsetrun=F)
  }
  
  #deal with the overlapping mss/tm composites
  print("dealing with any temporally overlapping MSS/TM composites")
  
  
  msscompfiles = list.files(mssdir, ".bsq$", recursive=T, full.names=T)
  tmcompfiles = list.files(tmdir, ".bsq$", recursive=T, full.names=T)
  thesetm = which(basename(tmcompfiles) %in% basename(msscompfiles))
  ref_files = sort(tmcompfiles[thesetm])
  thesemss = which(basename(msscompfiles) %in% basename(tmcompfiles))
  dep_files = sort(msscompfiles[thesemss])
  if(length(dep_files) > 0){
    for(i in 1:length(dep_files)){
      mssr= raster(dep_files[i])
      mssrnas = which(values(mssr)==0)
      mssr[mssrnas] = NA
      tmr= raster(ref_files[i])
      tmrnas = which(values(tmr)==0)
      tmr[tmrnas] = NA
      newimg = mosaic(mssr,tmr, fun="mean", na.rm=T)
      combnas = c(mssrnas, tmrnas)
      newimg[combnas] = 0
      projection(newimg) = set_projection(files[1])
      
      outimgfile = file.path(outdir,basename(dep_files[i]))
      writeRaster(newimg, outimgfile, format="ENVI", datatype = "INT2S",overwrite=T)
      envifilename = sub("bsq","envi",outimgfile)
      envixmlfile = paste(envifilename,".aux.xml",sep="")
      bsqxmlfile = sub("envi","bsq",envixmlfile)
      file.rename(envifilename,outimgfile)
      file.rename(envixmlfile,bsqxmlfile)
    }
  }
  
  print("dealing with any temporally overlapping ETM+/OLI composites")
  olicompfiles = list.files(olidir, ".bsq$", recursive=T, full.names=T)
  thesetm = which(basename(tmcompfiles) %in% basename(olicompfiles))
  ref_files = sort(tmcompfiles[thesetm])
  theseoli = which(basename(olicompfiles) %in% basename(tmcompfiles))
  dep_files = sort(olicompfiles[theseoli])
  if(length(dep_files) > 0){
    for(i in 1:length(dep_files)){
      mssr= raster(dep_files[i])
      mssrnas = which(values(mssr)==0)
      mssr[mssrnas] = NA
      tmr= raster(ref_files[i])
      tmrnas = which(values(tmr)==0)
      tmr[tmrnas] = NA
      newimg = mosaic(mssr,tmr, fun="mean", na.rm=T)
      combnas = c(mssrnas, tmrnas)
      newimg[combnas] = 0
      projection(newimg) = set_projection(files[1])
      
      outimgfile = file.path(outdir,basename(dep_files[i]))
      writeRaster(newimg, outimgfile, format="ENVI", datatype = "INT2S",overwrite=T)
      envifilename = sub("bsq","envi",outimgfile)
      envixmlfile = paste(envifilename,".aux.xml",sep="")
      bsqxmlfile = sub("envi","bsq",envixmlfile)
      file.rename(envifilename,outimgfile)
      file.rename(envixmlfile,bsqxmlfile)
    }
  }
    
  #rename files
  print("directory and file organization/cleaning")
  imglists = list.files(outdir, paste(runname,"_",index,"_composite_img_list.csv",sep=""), recursive=T, full.names=T)
  imglistyears = substr(basename(imglists),1,4)
  uniimglistyears = unique(imglistyears)
  for(i in 1:length(uniimglistyears)){
    outname = file.path(outdir,paste(uniimglistyears[i],"_",runname,"_",index,"_composite_img_list.csv", sep=""))
    theseones = which(imglistyears %in% uniimglistyears[i])
    if(length(theseones) == 1){file.rename(imglists[i],outname)}
    if(length(theseones) == 2){
      data1 = read.csv(imglists[theseones[1]])
      data2 = read.csv(imglists[theseones[2]])
      mergedlists = as.data.frame(rbind(data1,data2)$File)
      colnames(mergedlists) = "File"
      write.csv(mergedlists,outname)
    }
  }
  
  msstmolifiles = c(msscompfiles,tmcompfiles,olicompfiles)
  finalfiles = file.path(outdir,basename(msstmolifiles))
  
  for(i in 1:length(finalfiles)){
    check = file.exists(finalfiles[i])
    if(check == F){
      year = substr(basename(finalfiles[i]),1,4)
      files = list.files(dirname(msstmolifiles[i]),year,full.names=T)
      file.rename(files,file.path(outdir,basename(files)))
    }
  }
  
  #clean up
  unlink(c(mssdir,tmdir,olidir), recursive=T)
  
  
  print("making final annual composite stack")
  bname = paste(runname,"_",index,"_composite_stack.bsq", sep="")
  bands = sort(list.files(outdir, "composite.bsq$", full.names=T))
  fullnametif = file.path(outdir,bname)
  fullnamevrt = sub(".bsq", ".vrt", fullnametif)
  gdalbuildvrt(gdalfile=bands, output.vrt = fullnamevrt, separate=T) #, tr=c(reso,reso)
  gdal_translate(src_dataset=fullnamevrt, dst_dataset=fullnametif, of = "ENVI") #, co="INTERLEAVE=BAND"
  unlink(fullnamevrt)
}



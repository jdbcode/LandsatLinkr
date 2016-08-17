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
#' @import plyr
#' @export


mixel = function(msswrs1dir,msswrs2dir,tmwrs2dir,oliwrs2dir,index,outdir,runname,useareafile,doyears="all",order="none",overlap="mean",startday,endday,yearadj=0){
  
  mixel_find = function(files, refimg){
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
  
  mixel_mask = function(imgfile, useareafile, index){ #search,
    print(paste("...cloud masking:",basename(imgfile)))
    if(index == "tca" | index == "tcb"){band=1}
    if(index == "tcg"){band=2}
    if(index == "tcw"){band=3}
    
    sensor = substr(basename(imgfile), 1,2)
    if(sensor == "LM"){maskbit = "_cloudmask_30m.tif"} else {maskbit = "_cloudmask.tif"}
    
    maskfile = file.path(dirname(imgfile), paste(substr(basename(imgfile),1,16),maskbit,sep=""))
    img = raster(imgfile, band=band)
    mask = raster(maskfile)
    NAvalue(mask) = 0 #make 0 in the mask
    refimg = raster(useareafile)
    
    imgex = alignExtent(img, refimg, snap="near")
    maskex = alignExtent(mask, refimg, snap="near")
    extent(img) = imgex
    extent(mask) = maskex
    
    overlap = intersect(img, mask)
    mask = crop(mask, overlap)
    img = crop(img, overlap)
    
    img = img*mask
    return(img)
  }
  
  change_envi_to_bsq = function(file){
    envifilename = sub("bsq","envi",file)
    envixmlfile = paste(envifilename,".aux.xml",sep="")
    bsqxmlfile = sub("envi","bsq",envixmlfile)
    file.rename(envifilename,file)
    file.rename(envixmlfile,bsqxmlfile)
  }
  
  mixel_composite = function(outdir, imginfosub, runname, index, order, useareafile, overlap, yearadj){

    #outdir= mssdir 
    #imginfosub = mssdf

    uniyears = sort(unique(imginfosub$compyear))
    
    #for all the unique year make a composite
    for(i in 1:length(uniyears)){
      print(paste("working on year:", uniyears[i]))
      these = which(imginfosub$compyear == uniyears[i])
      theseimgs = imginfosub$file[these]

      if(order == "none"){theseimgs = theseimgs}

      len = length(theseimgs)
      #mask all the images in a year
      for(m in 1:len){
        mergeit = ifelse(m == 1, "r1", paste(mergeit,",r",m, sep=""))
        dothis = paste("r",m,"=mixel_mask(theseimgs[",m,"], useareafile, index)", sep="") 
        eval(parse(text=dothis))
      }
      
      #select a mosaic method
      if(overlap == "order"){mergeit = paste("newimg = merge(",mergeit,")", sep="")} else
      if(overlap == "mean"){mergeit = paste("newimg = mosaic(",mergeit,",fun=mean,na.rm=T)", sep="")} else
      if(overlap == "median"){mergeit = paste("newimg = mosaic(",mergeit,",fun=median,na.rm=T)", sep="")} else
      if(overlap == "max"){mergeit = paste("newimg = mosaic(",mergeit,",fun=max,na.rm=T)", sep="")} else
      if(overlap == "min"){mergeit = paste("newimg = mosaic(",mergeit,",fun=min,na.rm=T)", sep="")}
      
      #run the merge function
      print(paste("...merging files using ", overlap, ":",sep=""))
      for(h in 1:len){print(paste("......",basename(theseimgs[h]),sep=""))}
      if(len == 1){newimg = r1} else {eval(parse(text=mergeit))} #only run merge it if there are multiple files to merge
      
      #name the new file
      yearlabel = as.character(as.numeric(uniyears[i])+yearadj)
      newbase = paste(yearlabel,"_",runname,"_",index,"_composite.bsq", sep="")
      outimgfile = file.path(outdir,newbase)
      outtxtfile = sub("composite.bsq", "composite_img_list.csv", outimgfile)
      theseimgs = data.frame(theseimgs)
      colnames(theseimgs) = "File"
      write.csv(theseimgs, file=outtxtfile, row.names = F)
      
      #load in the usearea file and crop/extend the new image to it
      refimg = raster(useareafile)
      newimg = round(crop(newimg, refimg))
      newimg = extend(newimg, refimg, value=NA)

      #set NA values to 0 for use
      refimg = refimg != 0 # set all values not equal to 0 to 1 and 0 to 0 - NA can still be in there, it will be set to 0 in the final img
      refimg = refimg * newimg #set all 0's in the usearea file to 0 in the img file - if NA are in the usearea file, they will be transfered to img and then set to 0 in the final img
      newimg[is.na(newimg)] = 0 #set all NA to 0

      #write out the new image
      projection(newimg) = set_projection(files[1])
      writeRaster(newimg, outimgfile, format="ENVI", datatype = "INT2S",overwrite=T)
      change_envi_to_bsq(outimgfile)
    
      #clean the temp directory      
      delete_temp_files()
    }
  }
  
  delete_temp_files = function(){
    tempdir = dirname(rasterTmpFile())
    tempfiles = list.files(tempdir,full.names=T)
    unlink(tempfiles)
  }
  
  find_files = function(dir, search){
    if(length(which(is.na(dir) == T)) > 0){return(vector())} else{
      imgdir = normalizePath(file.path(dir,"images"),winslash="/")
      files = vector()
      for(i in 1:length(imgdir)){
        print(paste("finding files in: ",imgdir))
        files = c(files,list.files(imgdir[i], search, recursive=T, full.names=T))
      }
      if(length(files)==0){stop(
          paste("There were no tasselled cap files found in this directory: ",imgdir[i],".
          Make sure that you provided the correct scene head directory path and that all processing
          steps up to compositing have been completed. MSS directories should contain
          files with the extension 'tc_30m.tif' and 'tca_30m.tif', and TM/ETM+ and OLI 'tc.tif' and 'tca.tif'.
          A valid scene head directory path should look like this mock example: 'C:/mock/landsat/wrs2/045030'.
          The program will then append 'images' to the path and search recursively in that directory. If
          you wish to continue without this directory, re-run the compositing call and don't add this directory.",sep="")
        )
      } else{
        return(files)
      }
    }
  }
  
  combine_overlapping_senors = function(ref_files, dep_files){
    thesetm = which(basename(ref_files) %in% basename(dep_files))
    ref_files_sort = sort(ref_files[thesetm])
    thesemss = which(basename(dep_files) %in% basename(ref_files))
    dep_files_sort = sort(dep_files[thesemss])
    len = length(dep_files_sort)
    if(len > 0){
      for(i in 1:len){
        print(paste("...",i,"/",len,sep=""))
        depimg= raster(dep_files_sort[i])
        refimg= raster(ref_files_sort[i])
        NAvalue(depimg) = 0
        NAvalue(refimg) = 0
        newimg = mosaic(depimg,refimg, fun="mean", na.rm=T)
        newimg[is.na(newimg)] = 0
        projection(newimg) = set_projection(files[1])
        
        outimgfile = file.path(outdir,basename(dep_files_sort[i]))
        writeRaster(newimg, outimgfile, format="ENVI", datatype = "INT2S",overwrite=T)
        change_envi_to_bsq(outimgfile)
        
        delete_temp_files()
      }
    }
  }
  
  
  pixel_level_offset = function(ref_files, dep_files, outdir, sensor, projfile, runname){
    
    #find overlapping years
    theseref = which(basename(ref_files) %in% basename(dep_files))
    thesedep = which(basename(dep_files) %in% basename(ref_files))
    lendep = length(thesedep)
    lenref = length(theseref)
    if(lenref == 0 | lendep == 0 | lenref != lendep){return()} # get out if there are no matching files - no overlap
    
    print(paste("calculating pixel-level offset for: ",sensor," composites",sep=""))
    
    #sort the files to make sure they are in the same order
    ref_files_sort = sort(ref_files[theseref])
    dep_files_sort = sort(dep_files[thesedep])
    
    #find the mean pixel-wise difference between dependent and reference images
    for(i in 1:lendep){
      print(paste("...",basename(dep_files_sort[i]),sep=""))
      refimg = raster(ref_files_sort[i])
      depimg = raster(dep_files_sort[i])
      NAvalue(refimg) = 0
      NAvalue(depimg) = 0
      
      dif = refimg - depimg #get the difference 
      denom = !is.na(dif) #get the cells that are not NA after difference
      dif[is.na(dif)] = 0 #set the difference NA values to 0
      
      #make the sum difference and denominator layers
      if(i == 1){ #if i is one then start the layers
        difsum = dif
        denomsum = denom
      } else{ #else sum the layers
        difsum = sum(difsum, dif, na.rm=T)
        denomsum = sum(denomsum, denom, na.rm=T)
      }
    }
    print("...calculating mean pixel-level offset")
    meandiforig = round(difsum/denomsum) #mean the mean for the time series
    meandiforig[is.na(meandiforig)] = 0 #make division by 0 set to 0 instead of NA - there will be no correction for these pixels

    #figure out parts of file names
    if(sensor == "mss"){offsetdir = file.path(outdir,"mss_offset");deprunname="lm"; refrunname="lt"; meandiffilebname = "mss_mean_dif.bsq"}
    if(sensor == "oli"){offsetdir = file.path(outdir,"oli_offset");deprunname="lc"; refrunname="le"; meandiffilebname = "oli_mean_dif.bsq"}
    
    #write out the mean pixel-wise difference file
    dir.create(offsetdir, recursive=T, showWarnings=F)
    projection(meandiforig) = set_projection(projfile)
    meandiffile = file.path(offsetdir,meandiffilebname)
    writeRaster(meandiforig, meandiffile, format="ENVI", datatype = "INT2S",overwrite=T)
    change_envi_to_bsq(meandiffile)
    
    #write out the frequency file
    projection(denomsum) = set_projection(projfile)
    denomsumfile = file.path(dirname(meandiffile),"overlap_frequency.bsq")
    writeRaster(denomsum, denomsumfile, format="ENVI", datatype="INT2S",overwrite=T)
    change_envi_to_bsq(denomsumfile)
    
    #move files around
    for(i in 1:lendep){
      from_dep_files = list.files(dirname(dep_files_sort[1]),substr(basename(dep_files_sort[i]),1,4), full.names = T)
      from_ref_files = list.files(dirname(ref_files_sort[1]),substr(basename(ref_files_sort[i]),1,4), full.names = T)
      to_dep_files = file.path(offsetdir,sub(paste("_",runname,"_",sep=""),paste("_",deprunname,"_",sep=""),basename(from_dep_files)))
      to_ref_files = file.path(offsetdir,sub(paste("_",runname,"_",sep=""),paste("_",refrunname,"_",sep=""),basename(from_ref_files)))
      
      file.copy(from_dep_files, to_dep_files)
      file.copy(from_ref_files, to_ref_files)
    }
    
    #adjust the dep images
    print("...adjusting images by mean pixel-level offset:")
    for(i in 1:length(dep_files)){
      print(paste("......",basename(dep_files[i]),sep=""))
      img = raster(dep_files[i])
      NAvalue(img) = 0
      
      meandiforig = raster(meandiffile) #need to load this each time because delete_temp_files() gets called at the end of each loop - if this file is big it is held in the temp directory and will be deleted
      img = img + meandiforig
      img[is.na(img)] = 0
      
      #write out the new image
      projection(img) = set_projection(projfile)
      
      #delete the old .bsq files - keep the "img_list.csv" files
      allfiles = list.files(dirname(dep_files[i]),substr(basename(dep_files[i]),1,nchar(basename(dep_files[i]))-4),full.names=T)
      allfiles = grep("img_list.csv", allfiles, invert=T, value=T) #keep the "img_list.csv" files
      unlink(allfiles)
      writeRaster(img, dep_files[i], format="ENVI", datatype="INT2S",overwrite=T)
      change_envi_to_bsq(dep_files[i])
      
      delete_temp_files()
    }
  }
  
  #check for leap year
  leapyear = function(year){
    return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
  }
  
  #create decimal year day
  decyearday = function(year,day){
    num = ifelse(leapyear(year) == T, 366, 365)
    return(year+day/num)
  }
  
  #create decimal day
  decday = function(day){
    num = ifelse(day == 366, 366, 365)
    return(day/num)
  }

  #########################################################################
  #########################################################################
  #########################################################################
  print(paste("working on index:",index))
  
  #create some search terms depending on index
  if(index == "tca"){msssearch="tca_30m.tif$"; tmsearch="tca.tif$"; olisearch="tca.tif$"}
  if(index == "tcb"){msssearch="tc_30m.tif$"; tmsearch="tc.tif$"; olisearch="tc.tif$"}
  if(index == "tcg"){msssearch="tc_30m.tif$"; tmsearch="tc.tif$"; olisearch="tc.tif$"}
  if(index == "tcw"){msssearch="tc_30m.tif$"; tmsearch="tc.tif$"; olisearch="tc.tif$"}
  
  #find the files
  msswrs1files = find_files(msswrs1dir, msssearch)
  msswrs2files = find_files(msswrs2dir, msssearch)
  tmwrs2files = find_files(tmwrs2dir, tmsearch)
  oliwrs2files = find_files(oliwrs2dir, olisearch)
  
  #put all the files together in a vector and check to make sure the files intersect the useareafile, if not they will be excluded
  files = c(msswrs1files,msswrs2files,tmwrs2files,oliwrs2files)
  
  #find files that intersect the usearea file
  files = mixel_find(files, useareafile)
  if(length(files)==0){stop("There were no files in the given directories that intersect the provided 'usearea file'.
                            Make sure that you provided the correct file, checked that it actually overlaps the scenes you 
                            specified for compositing, and that the projection is the same as the images.")}
  
  
  #create a table with info on the files
  imginfo = data.frame(file = as.character(files))
  imginfo$file = as.character(imginfo$file)
  bname = basename(imginfo$file)
  imginfo$year = substr(bname, 10, 13)
  imginfo$day = substr(bname, 14,16)
  imginfo$sensor = substr(bname, 1,2)
  imginfo$compyear = imginfo$decdate = NA
  for(i in 1:nrow(imginfo)){imginfo$decdate[i] = decyearday(as.numeric(imginfo$year[i]),as.numeric(imginfo$day[i]))}
  
  #figure out the dec year day range that is good
  #uniyears = as.numeric(sort(unique(imginfo$year)))
  uniyears = 1972:as.numeric(format(Sys.Date(),'%Y'))
  
  if(doyears != "all"){uniyears = uniyears[match(doyears,uniyears)]}

  decstart = decday(startday)
  decend = decday(endday)
  dif = ifelse(decstart > decend, (1-decstart)+decend, decend - decstart)
  start = uniyears+decstart
  end = start+dif
  yearsdf = data.frame(uniyears,start,end)
  
  #figure out which images are in the composite date range, get rif of ones that aren't
  for(i in 1:nrow(yearsdf)){
    these = which(imginfo$decdate >= yearsdf$start[i] & imginfo$decdate <= yearsdf$end[i])
    imginfo$compyear[these] = yearsdf$uniyears[i]
  }
  imginfo = na.omit(imginfo)
  if(nrow(imginfo)==0){stop("There were no files in the given directories that intersect the provided start and end year-of-day bounds.
                            Make sure that you provided the correct limits and check that there are actually image files that intersect 
                            the specified date range.")}
  
  
  #make composites
  mssdir = file.path(outdir,"mss")
  tmdir = file.path(outdir,"tm")
  olidir = file.path(outdir,"oli")
  
  #pull out files by sensor
  mssdf = imginfo[imginfo$sensor == "LM",]
  tmetmdf = imginfo[imginfo$sensor == "LT" | imginfo$sensor == "LE",]
  olidf = imginfo[imginfo$sensor == "LC",]
  
  #create annual composites for all sensors
  if(nrow(mssdf) != 0){
    dir.create(mssdir, recursive=T, showWarnings=F)
    mixel_composite(mssdir, mssdf, runname=runname,index=index, order=order, useareafile=useareafile, overlap=overlap, yearadj=yearadj)
  }
  if(nrow(tmetmdf) != 0){
    dir.create(tmdir, recursive=T, showWarnings=F)
    mixel_composite(tmdir, tmetmdf, runname=runname,index=index, order=order, useareafile=useareafile, overlap=overlap, yearadj=yearadj)
  }
  if(nrow(olidf) != 0){
    dir.create(olidir, recursive=T, showWarnings=F)
    mixel_composite(olidir, olidf, runname=runname,index=index, order=order, useareafile=useareafile, overlap=overlap, yearadj=yearadj)
  }
  
  #deal with the overlapping composites
  msscompfiles = list.files(mssdir, ".bsq$", recursive=T, full.names=T)
  tmcompfiles = list.files(tmdir, ".bsq$", recursive=T, full.names=T)
  olicompfiles = list.files(olidir, ".bsq$", recursive=T, full.names=T)
  
  pixel_level_offset(ref_files=tmcompfiles, dep_files=msscompfiles, outdir=outdir, sensor="mss", projfile=files[1], runname=runname)
  pixel_level_offset(ref_files=tmcompfiles, dep_files=olicompfiles, outdir=outdir, sensor="oli", projfile=files[1], runname=runname)
  
  
  print("dealing with any temporally overlapping MSS/TM composites")
  combine_overlapping_senors(tmcompfiles, msscompfiles)
  
  print("dealing with any temporally overlapping ETM+/OLI composites")
  combine_overlapping_senors(tmcompfiles, olicompfiles)
  
  
  #rename files
  print("directory and file organization/cleaning")
  imglists = list.files(outdir, paste(runname,"_",index,"_composite_img_list.csv",sep=""), recursive=T, full.names=T)
  imglistyears = substr(basename(imglists),1,4)
  uniimglistyears = unique(imglistyears)
  for(i in 1:length(uniimglistyears)){
    outname = file.path(outdir,paste(uniimglistyears[i],"_",runname,"_",index,"_composite_img_list.csv", sep=""))
    theseones = which(imglistyears %in% uniimglistyears[i])
    for(l in 1:length(theseones)){
      if(l == 1){
        alllimglist = read.csv(imglists[theseones[l]], stringsAsFactors=F)$File
      } else{
        alllimglist = c(alllimglist, read.csv(imglists[theseones[l]], stringsAsFactors=F)$File)
      }
    }
    alllimglist = data.frame(File = alllimglist)
    write.csv(alllimglist, outname, row.names = F)
  }
  
  #move files
  msstmolifiles = c(msscompfiles, tmcompfiles, olicompfiles)
  finalfiles = file.path(outdir, basename(msstmolifiles))
  for(i in 1:length(finalfiles)){
    check = file.exists(finalfiles[i])
    if(check == F){
      year = substr(basename(finalfiles[i]), 1, 4)
      files = list.files(dirname(msstmolifiles[i]), year, full.names=T)
      file.rename(files, file.path(outdir, basename(files)))
    }
  }
  
  #clean up
  unlink(c(mssdir,tmdir,olidir), recursive=T)
  
  #make the final stack
  print("making final annual composite stack")
  bname = paste(runname,"_",index,"_composite_stack.bsq", sep="")
  bands = sort(list.files(outdir, "composite.bsq$", full.names=T))
  fullnametif = file.path(outdir,bname)
  fullnamevrt = change_extension("bsq", "vrt", fullnametif)
  gdalbuildvrt(gdalfile=bands, output.vrt = fullnamevrt, separate=T) #, tr=c(reso,reso)
  gdal_translate(src_dataset=fullnamevrt, dst_dataset=fullnametif, of = "ENVI") #, co="INTERLEAVE=BAND"
  unlink(fullnamevrt)
}



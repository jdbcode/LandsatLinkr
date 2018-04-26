#' Spatially warp an MSS image
#'
#' Spatially warp an MSS image to match the spatial properties of a reference image
#' @param reffile character. MSS image file that has low spatial RMSE and low cloud cover
#' @param fixfile character. MSS image file to be spatially warped to match the reference file
#' @param window numeric. image subset size used to define cross-correlation calculation. unit is pixels along one side of a square
#' @param search numeric. neighborhood search window size in which to find tie-point offset. unit is pixels along one side of a square
#' @param sample numeric. target number of tie-points
#' @param refstart numeric. c(xcoord,ycoord). reference image coordinate for a pixel identified as common in both the reference and the to-be-warped image. used to calculate an initial offset between the two images.
#' @param fixstart numeric. c(xcoord,ycoord). fix image coordinate for a pixel identified as common in both the reference and the to-be-warped image. used to calculate an initial offset between the two images.
#' @import raster
#' @import gdalUtils
#' @export


msswarp = function(reffile, fixfile, refstart=c(0,0), fixstart=c(0,0)){

  mode = 'warp'
  method = 'order 2'
  
  #mode can be: "rmse" or "warp"
  #method can be: "tps" or "order 1" or "order 2"
  
  
  #set default parameters
  search=35 #27 
  if(mode == "rmse"){search = 27} #it should be pretty close, no need to look over a larger region - keep at 27 since that is what the thresholds were selected at.
  
  #scale image values to center on mean and 1 unit variance (global mean and variance)
  scaleit = function(matrix){
    stnrd = (matrix - (mean(matrix, na.rm = TRUE)))/(sd(matrix, na.rm = TRUE))
    return(stnrd)
  }
  
  #make a kernal around a given point
  make_kernal = function(img, point1, windowsize){
    radius = floor(windowsize/2)
    ccol = colFromX(img, point1[1])
    crow = rowFromY(img, point1[2])
    mincol = ccol-radius
    maxcol = ccol+radius
    minrow = crow-radius
    maxrow = crow+radius
    return(extent(c(mincol,maxcol,minrow,maxrow)))
  }
  
  calc_rmse = function(info,reso){
    xresid = (info[,"refx"]-info[,"fixx"])^2 #get the residuals of each x
    yresid = (info[,"refy"]-info[,"fixy"])^2 #get the residuals of each y
    r = (sqrt(xresid+yresid))/reso #get the rmse of each xy point
    x_rmse = sqrt(mean(xresid))/reso
    y_rmse = sqrt(mean(yresid))/reso
    total_rmse = sqrt((x_rmse^2)+(y_rmse^2)) #total rmse including all points
    rmse_info = list(x_rmse=x_rmse, y_rmse=y_rmse, total_rmse=total_rmse, r=r)
    return(rmse_info)
  }
  
  make_gdaltrans_opts = function(info, wktfile, fixfile, tempname){
    info[,"refx"] = info[,"refx"]+(reso/2)
    info[,"refy"] = info[,"refy"]-(reso/2)
    fixcol = paste(info[,"fixcol"]) #fix col for tie point
    fixrow = paste(info[,"fixrow"]) #fix row for tie point
    refx = paste(round(info[,"refx"]))  #fix x tie point coord 
    refy = paste(round(info[,"refy"]))  #fix y tie point coord
    gcpstr = paste(" -gcp", fixcol, fixrow, refx, refy, collapse="")
    gdaltrans_cmd = paste("-of Gtiff -ot Byte -co INTERLEAVE=BAND -a_srs", wktfile, gcpstr, fixfile, tempname)
    return(gdaltrans_cmd)
  }
  
  run_the_file = function(fixfile){
    info = get_metadata(fixfile)
    dt = as.character(info$datatype)
    rmsefile = sub("archv.tif","cloud_rmse.csv",fixfile)
    rmse = as.numeric(read.table(rmsefile,sep=",")[3])
    runit = as.numeric(rmse > 0.75 & dt == "L1T")
    return(runit)
  }
  
  
  
  runit = run_the_file(fixfile)
  if(runit == 0){return(0)}
  
  print(paste('Working on image:', basename(fixfile)))
  
  #read in the fix image
  #fiximg = raster(fixfile, band=4) #load the fix image
  fiximg = brick(fixfile) #load the fix image
  origfiximg = subset(fiximg, subset=4) #save a copy of the original fix image so that the tie point coords can be assigned to the original row and cols
  
  #shift the fiximg if there is an initial offset provided
  shiftit = refstart - fixstart
  if(sum(shiftit) != 0){fiximg = shift(fiximg, x=shiftit[1], y=shiftit[2])}
  
  #load the ref image
  refimg = raster(reffile, 4) 
  
  #make sure that the ref and fix img are croppd to eachother
  extent(fiximg)  = alignExtent(fiximg, refimg, snap="near")
  ext = intersect(extent(refimg), extent(fiximg))
  refimg = crop(refimg, ext)
  fiximg = crop(fiximg, ext)
  #refimg = intersect(refimg, fiximg)
  #fiximg = intersect(fiximg, refimg)
  
  #fiximgb1 = raster(fixfile, band=1)
  #fiximgb4 = raster(fixfile, band=4)
  #get bands 1 and 4 out for cloud and shadow id in the fix image, as well as for cross correlation with the reference image
  fiximgb1 = subset(fiximg, subset=1)
  fiximgb4 = subset(fiximg, subset=4)
  fiximg = fiximgb4 #need to copy because fiximg will be scaled, but we also need an unalted copy to find shadows in
  
  
  #calculate similarity index input values from the fix image subset
  values(fiximg) = scaleit(values(fiximg))
  values(fiximg)[is.na(values(fiximg))] = 0
  
  #calculate similarity index input values from the ref image subset
  values(refimg) = scaleit(values(refimg))
  values(refimg)[is.na(values(refimg))] = 0
  refimgsqr = refimg^2
  
  #get the resolution  
  reso = xres(refimg)
  
  #adjust the window and search size so that they are odd numbers
  if (search %% 2 == 0){search = search+1}
  
  #sample the reference image, laying down a regular grid of points to check
  s = sampleRegular(refimg, 15000, cells=T)[,1]
  xy = xyFromCell(refimg,s) #[,1] #get the xy coordinates for each good point
  
  #filter points in fiximg that fall on clouds
  thesecells = na.omit(cellFromXY(fiximgb1, xy))   #get fiximg cell index for sample 
  b = which(fiximgb1[thesecells] < 120 & 
            fiximgb4[thesecells] > 30 &
            refimg[thesecells] > 0) # fiximgb1[theseones] != NA) #exclude points that don't meet criteria
  
  #if the number of sample points is less than 10 delete the image return
  
  #TODO - get out if there are not enough points to work with - there are examples of for doing this below
  #if(length(b) < 10){
  #  delete_files(fixfile, 2)
  #  return(0)
  #}
  
  #subset the sample
  n_samp = length(b)
  cat(paste("n points from original sample:",n_samp),"\n")
  n_subsamp = 1200
  if(n_samp < n_subsamp){n_subsamp = n_samp}
  subsamp = sample(b, n_subsamp)
  xy = xy[subsamp,] #subset the original sample
  s = s[subsamp] #subset the original sample
  rowcol = rowColFromCell(refimg, s)
  
  #make an array to hold all information collected
  info = cbind(c(0),xy, rowcol[,2], rowcol[,1], array(0, c(length(xy[,1]), 8)))
  cnames = c("point","refx","refy","refcol","refrow","fixx","fixy","fixcol","fixrow","nmax", "max","edgedist","decision")
  colnames(info) = cnames
  
  #iterate process of creating a similarity surface for each check point in the sample
  window_size = c(101,201,275)
  for(size in 1:3){
    cat(paste("Working on window size set: ",size,"/3",sep=""),"\n")
    if(mode != "rmse"){
      if(size == 1){pdf_file = sub("archv.tif", "ccc_surface_100w.pdf",fixfile)}
      if(size == 2){pdf_file = sub("archv.tif", "ccc_surface_200w.pdf",fixfile)}
      if(size == 3){pdf_file = sub("archv.tif", "ccc_surface_275w.pdf",fixfile)}
      unlink(pdf_file) #delete the pdf if it exists
      pdf(file=pdf_file, width=10, heigh=7) #size of the pdf page
      par(mfrow=c(2,3)) #number of trajectories to place on a page (columns, rows)
    }
    window = window_size[size]
    
    #adjust the window and search size so that they are odd numbers
    if (window %% 2 == 0){window = window+1}
    radius = floor(window/2) #radius of the window in pixels
    nrc = search+(radius*2) #the reference extent length to slide over
    
    
    for(point in 1:length(info[,1])){ 
      #print(point) #print the point so we know where we're at
      if(info[point,"decision"] == 1){
        #  print("already good, skipping...")
        next
      }
      if(size == 1){info[point,"point"] = point} #info[point,1] = point #put the point number into the info table
      
      #make a subset of the reference image for the fiximg chip to slide over
      a = make_kernal(refimg, info[point,2:3], nrc)
      test = c(a@ymax,a@ymin,a@xmin,a@xmax)
      
      if(sum(is.na(test)) > 0){next}
      if(sum(test < 0) > 0){next}
      if(a@ymax > nrow(refimg) | a@ymin > nrow(refimg)){next}
      if(a@xmax > ncol(refimg) | a@xmin > ncol(refimg)){next}
      ext=extent(refimg,a@ymin,a@ymax,a@xmin,a@xmax)
      refsub = crop(refimg, ext)
      
      #make subset of fiximg (fiximg chip)
      a = make_kernal(fiximg, info[point,2:3], window)
      test = c(a@ymax,a@ymin,a@xmin,a@xmax)
      if(sum(is.na(test)) > 0){next}
      if(sum(test < 0) > 0){next}
      if(a@ymax > nrow(fiximg) | a@ymin > nrow(fiximg)){next}
      if(a@xmax > ncol(fiximg) | a@xmin > ncol(fiximg)){next}
      ext=extent(fiximg,a@ymin,a@ymax,a@xmin,a@xmax)
      fixsub = crop(fiximg, ext)
      
      #create numerator
      tofix = matrix(values(fixsub),ncol=window,byrow = T)
      
      if(length(tofix) %% 2 == 0) {
        #cat("Skipping","\n")
        next
      }
      
      num = focal(refsub, w=tofix ,fun=sum)
      
      #get refimg denom
      a = make_kernal(refimgsqr, info[point,2:3], nrc)
      ext=extent(refimgsqr,a@ymin,a@ymax,a@xmin,a@xmax)
      refsubsqr = crop(refimgsqr, ext)
      sumrefsubsqr = focal(refsubsqr, w=matrix(1,window, window)) #get the summed product of the refsubimg
      sumfixsubsqr = sum(values(fixsub)^2) #fiximg standard only gets calcuated once
      denom = sqrt(sumfixsubsqr*sumrefsubsqr)
      
      #badone=0
      if(cellStats(num, stat="sum") + cellStats(denom, stat="sum") == 0){next} 
      
      ncc = num/denom
      buf = (nrow(ncc)-search)/2
      off1 = buf+1
      off2 = buf+search
      ext = extent(ncc,off1,off2,off1,off2)
      ncc = crop(ncc,ext)
      nccv = values(ncc)
      nccm = matrix(nccv, ncol=sqrt(length(nccv)), byrow=T)
      
      x = y = seq(1,ncol(nccm),1)
      
      good = which(values(ncc) == maxValue(ncc))[1] 
      
      ####
      xoffsetcoord = xFromCell(ncc, good)
      yoffsetcoord = yFromCell(ncc, good)
      xoffset = xoffsetcoord - info[point,"refx"]
      yoffset = yoffsetcoord - info[point,"refy"]
      info[point,"fixx"] = xoffsetcoord-(xoffset*2)
      info[point,"fixy"] = yoffsetcoord-(yoffset*2)
      ####
      
      #get the row and column numbers for the fix image
      origfiximg_x = info[point,"fixx"]-shiftit[1]
      origfiximg_y = info[point,"fixy"]-shiftit[2]
      a = cellFromXY(origfiximg, c(origfiximg_x,origfiximg_y))
      fiximgrc = rowColFromCell(origfiximg, a)
      info[point,"fixcol"] = fiximgrc[2] #info[point,8] = fiximgrc[2]
      info[point,"fixrow"] = fiximgrc[1] #info[point,9] = fiximgrc[1]
      
      
      #############################################################
      #screen by outlier
      ccc = scaleit(nccm)  
      maxmat = max(ccc, na.rm = T)
      rowcol = which(ccc == maxmat, arr.ind=T)
      r = ccc[rowcol[1],] 
      c = ccc[,rowcol[2]] 
      rmax = which(r == maxmat) #max(r)
      cmax = which(c == maxmat) #max(c)
      
      dist1 = abs(c((rmax - 1), (rmax - length(r)), (cmax - 1), (cmax - length(c))))/(floor(search/2)) 
      dist2 = min(dist1)
      
      #place filtering values in info table
      info[point,"nmax"] = length(which(ccc == maxmat)) #2 #info[point,10] = length(which(ccc == maxmat)) #2
      info[point,"max"] = round(maxmat, digits=1) #3 #info[point,11] = round(maxmat, digits=1) #3
      info[point,"edgedist"] = round(dist2, digits=2) #6 #info[point,14] = round(dist2, digits=2) #6
      
      #decide what surfaces are good\bad
      bad = array(0,3)
      bad[1] = info[point,"nmax"] > 1 #number of max peaks eq 1
      bad[2] = info[point,"max"] < 3 #peak ge to 3 standard devs from mean
      bad[3] = info[point,"edgedist"] < 0.12 #peak distance from edge >= 0.12
      info[point,"decision"] = sum(bad) == 0 #7
      #if(badone == 1){info[point,"decision"] = 0}
      
      #filter plots that will crash the plotting because of weird data points (na, NaN, (-)Inf)
      bad1 = is.na(ccc)
      bad2 = is.infinite(ccc) 
      bad3 = is.nan(ccc)
      badindex = which(bad1 == T | bad2 == T | bad3 == T)
      ccc[badindex] = 0
      if(length(which(ccc == 0)) == length(ccc)){next}
      
      #plot the cross correlation surface
      #title = paste(point,info$nmax[point],info$max[point],info$edgedist[point], info$decision[point], sep = ",")
      if(info[point,"decision"] == 1){status = "accept"} else {status = "reject"}
      title = paste("Point:", point, status)
      #print(title)
      if(mode != "rmse"){
        ccc = ccc[nrow(ccc):1,]
        persp(x, y, ccc, theta = 30, phi = 30, expand = 0.5, col = 8, main=title)
      }
    }
    cat(paste("n goods =",length(which(info[,"decision"] == 1))),"\n")
    if(mode != "rmse"){
      dev.off() #turn off the plotting device
    }
  }
  #print(paste("n goods =",length(which(info[,"decision"] == 1))))
  
  #write all the point info to a file
  if(mode != "rmse"){
    info_file = sub("archv.tif", "info_full.csv",fixfile)
    write.csv(info, file=info_file, row.names = F) 
  }
  
  #subset the points that passed the surface tests
  these = which(info[,"decision"] == 1)
  
  # make file info for rmse - not used if not an rmse run
  rmse_outfile = file.path(dirname(fixfile),paste(substr(basename(fixfile),1,16),"_rmse.Rdata",sep=""))
  rmse_info = list(calc_rmse = F, x_rmse=NA, y_rmse=NA, total_rmse=NA, info=info)  
  
  
  
  #if the number of sample points is less than 10 delete the image and return
  if(length(these) < 10){
    delete_files(fixfile, 2)
    return(0)
  }
  
  info = info[these,]
  
  #filter points based on rmse contribution
  if(mode != "rmse"){
    rmse = calc_rmse(info,reso)
    r = rmse$r
    sdr = sd(r)
    meanr = mean(r)
    limit = meanr+sdr*2
    goods = which(r <= limit)
    n_outliers = nrow(info)-length(goods)
    info = info[goods,]
    cat(paste("Getting rid of:",n_outliers,"outliers"),"\n")
    cat(paste("There are still:",nrow(info),"points"),"\n")
    #maxr = endit = 10
    #while(maxr >2 & endit != 0){
    #  rmse = calc_rmse(info,reso)
    #  if (rmse$total_rmse != 0){contr = rmse$r/rmse$total_rmse} else contr = rmse$r #error contribution of each point
    #  maxr = max(contr) #while loop controler
    #  b = which(contr < 2) #subset finder - is point 2 times or greater in contribution
    #  info = info[b,] #subset the info based on good rsme
    #  endit = sum(contr[b]) #while loop controler
    #}
  }
  
  #if the number of sample points is less than 10 delete the image and return
  if(length(these) < 10){
    delete_files(fixfile, 2)
    return(0)
  }
  
  #if this is an rmse run, then save the info and get out
  if(mode == "rmse"){
    rmse = calc_rmse(info,reso)
    
    rmse_info$calc_rmse = T
    rmse_info$x_rmse=rmse$x_rmse
    rmse_info$y_rmse=rmse$y_rmse
    rmse_info$total_rmse=rmse$total_rmse
    rmse_info$info=info
    
    save(rmse_info, file=rmse_outfile)
    return(0)
  }
  
  #write out the filtered points that will be used in the transformation
  info_file = sub("archv.tif", "info_sub.csv",fixfile)
  write.csv(info, file=info_file, row.names = F)
  
  #make some output file names
  tempname = sub("archv.tif", "temp.tif", fixfile) #"K:/scenes/034032/images/1976/LM10360321976248_archv_l1g_warp.tif"
  #outfile = sub("archv.tif", "archv_l1g2l1t.tif", fixfile)
  wktfile = sub("archv.tif","wkt.txt", fixfile)
  #gcpfile = sub("archv.tif", "gcp.txt", fixfile)
  gdaltransoptsfile = sub("archv.tif", "gdal_trans_opts.txt", fixfile)
  
  #write out a projection file for gdal translate to use
  proj = system(paste("gdalsrsinfo -o wkt", fixfile), intern = TRUE)
  write(proj, wktfile)
  
  #create the warp cmd and save as a file
  gdaltrans_opts = make_gdaltrans_opts(info, wktfile, fixfile, tempname)
  write(gdaltrans_opts, file=gdaltransoptsfile)
  
  #run the warp command file
  cmd = paste("gdal_translate --optfile", gdaltransoptsfile)
  system(cmd)
  
  #gdal warp command
  gdalwarp_cmd = paste("gdalwarp -of Gtiff", paste("-", method, sep=""), "-ot Byte -srcnodata 0 -dstnodata 0 -co INTERLEAVE=BAND -overwrite -multi -tr", reso, reso, tempname, fixfile) #fixfile   "-tps"  "-order 2", "-order 3" 
  system(gdalwarp_cmd)
  
  
  # ######################## warping method tests#####################################################
  # outfiletest = sub("archv.tif", "archv_l1g2l1t_test_tps.tif", fixfile)
  # gdalwarp_cmd = paste("gdalwarp -of Gtiff -tps -ot Byte -srcnodata 0 -dstnodata 0 -co INTERLEAVE=BAND -overwrite -multi -tr", reso, reso, tempname, outfiletest) #fixfile   "-tps"  "-order 2", "-order 3" 
  # system(gdalwarp_cmd)
  # 
  # outfiletest = sub("archv.tif", "archv_l1g2l1t_test_order1.tif", fixfile)
  # gdalwarp_cmd = paste("gdalwarp -of Gtiff -order 1 -ot Byte -srcnodata 0 -dstnodata 0 -co INTERLEAVE=BAND -overwrite -multi -tr", reso, reso, tempname, outfiletest) #fixfile   "-tps"  "-order 2", "-order 3" 
  # system(gdalwarp_cmd)
  # 
  # outfiletest = sub("archv.tif", "archv_l1g2l1t_test_order2.tif", fixfile)
  # gdalwarp_cmd = paste("gdalwarp -of Gtiff -order 2 -ot Byte -srcnodata 0 -dstnodata 0 -co INTERLEAVE=BAND -overwrite -multi -tr", reso, reso, tempname, outfiletest) #fixfile   "-tps"  "-order 2", "-order 3" 
  # system(gdalwarp_cmd)
  # ##################################################################################################
  
  
  #delete the temp file
  unlink(list.files(dirname(fixfile), pattern = "temp", full.names = T))
  return(1)
}
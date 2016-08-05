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

msswarp = function(reffile, fixfile, window=275, search=27, sample=1000, refstart=c(0,0), fixstart=c(0,0)){
  
  
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
    minrow = crow+radius
    maxrow = crow-radius
    return(extent(c(mincol,maxcol,minrow,maxrow)))
  }
   
  #check to see if the image should be run
  info = get_metadata(fixfile)
  dt = as.character(info$datatype)
  rmsefile = sub("archv.tif","cloud_rmse.csv",fixfile)
  rmse = as.numeric(read.table(rmsefile,sep=",")[3])
  runit = as.numeric(rmse > 0.5 & dt == "L1T")
  if(runit == 1){
    #read in the fix image
    fiximg = raster(fixfile, band=3) #load the fix image
    origfiximg = fiximg #save a copy of the original fix image
    fiximgb1 = raster(fixfile, band=1)
    
    #shift the fiximg if there is an initial offset provided
    shiftit = refstart - fixstart
    if(sum(shiftit) != 0){fiximg = shift(fiximg, x=shiftit[1], y=shiftit[2])}
    
    #load the ref image
    refimg = raster(reffile, 3) 
    
    #make sure that the ref and fix img are croppd to eachother
    refimg = intersect(refimg, fiximg)
    fiximg = intersect(fiximg, refimg)
    
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
    if (window %% 2 == 0){window = window+1}
    if (search %% 2 == 0){search = search+1}
    radius = floor(window/2) #radius of the window in pixels
    nrc = search+(radius*2) #the reference extent length to slide over
    
    #sample the the reference image, laying down a regular grid of points to check
    s = sampleRegular(refimg, sample, cells=T)
    s = s[,1]
    xy = xyFromCell(refimg,s) #[,1] #get the xy coordinates for each good point
    
    #filter points in fiximg that fall on clouds
    theseones = cellFromXY(fiximgb1, xy)   #get fiximg cell index for sample 
    theseones = na.omit(theseones)
    a = fiximgb1[theseones] #extract values for fiximg cell sample
    b = which(fiximgb1[theseones] < 100) # fiximgb1[theseones] != NA) #exclude points that don't meet criteria
    
    #if the number of sample points is less than 10 delete the image return
    if(length(b) < 10){
      delete_files(fixfile, 2)
      return(0)
    }
    
    #subset the sample
    xy = xy[b,] #subset the original sample
    s = s[b] #subset the original sample
    rowcol = rowColFromCell(refimg, s)
    
    #make an array to hold all information collected
    info = cbind(c(0),xy, rowcol[,2], rowcol[,1], array(0, c(length(xy[,1]), 8)))
    cnames = c("point","refx","refy","refcol","refrow","fixx","fixy","fixcol","fixrow","nmax", "max","edgedist","decision")
    colnames(info) = cnames
    
    #start a pdf file to hold image chips and similarity surfaces
    pdf_file = sub("archv.tif", "ccc_surface.pdf",fixfile)
    unlink(pdf_file) #delete the pdf if it exists
    pdf(file=pdf_file, width=10, heigh=7) #size of the pdf page
    par(mfrow=c(2,3)) #number of trajectories to place on a page (columns, rows)

    #iterate process of creating a similarity surface for each check point in the sample
    for(point in 1:length(info[,1])){ 
      print(point) #print the point so we know where we're at
      info[point,"point"] = point #info[point,1] = point #put the point number into the info table
      
      #make a subset of the reference image for the fiximg chip to slide over
      a = make_kernal(refimg, info[point,2:3], nrc)
      test = c(a@ymax,a@ymin,a@xmin,a@xmax)
      
      if(sum(is.na(test)) > 0){next}
      if(sum(test < 0) > 0){next}
      if(a@ymax > nrow(refimg) | a@ymin > nrow(refimg)){next}
      if(a@xmax > ncol(refimg) | a@xmin > ncol(refimg)){next}
      ext=extent(refimg,a@ymax,a@ymin,a@xmin,a@xmax)
      refsub = crop(refimg, ext)
      
      #make subset of fiximg (fiximg chip)
      a = make_kernal(fiximg, info[point,2:3], window)
      test = c(a@ymax,a@ymin,a@xmin,a@xmax)
      if(sum(is.na(test)) > 0){next}
      if(sum(test < 0) > 0){next}
      if(a@ymax > nrow(fiximg) | a@ymin > nrow(fiximg)){next}
      if(a@xmax > ncol(fiximg) | a@xmin > ncol(fiximg)){next}
      ext=extent(fiximg,a@ymax,a@ymin,a@xmin,a@xmax)
      fixsub = crop(fiximg, ext)
      
      #create numerator
      tofix = matrix(values(fixsub),ncol=window,byrow = T)
      
      if (length(tofix) %% 2 == 0) {
        print("skipping")
        next
      }
      
      num = focal(refsub, w=tofix ,fun=sum)
      
      #get refimg denom
      a = make_kernal(refimgsqr, info[point,2:3], nrc)
      ext=extent(refimgsqr,a@ymax,a@ymin,a@xmin,a@xmax)
      refsubsqr = crop(refimgsqr, ext)
      sumrefsubsqr = focal(refsubsqr, w=matrix(1,window, window)) #get the summed product of the refsubimg
      sumfixsubsqr = sum(values(fixsub)^2) #fiximg standard only gets calcuated once
      denom = sqrt(sumfixsubsqr*sumrefsubsqr)
      
      badone=0
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
      
      good = which(values(ncc) == maxValue(ncc)) 
      good = good[1]
      
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
      #a = cellFromXY(origfiximg, c(info[point,"refx"],info[point,"refy"])) #fiximg
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
      if(badone == 1){info[point,"decision"] = 0}
      
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
      persp(x, y, ccc, theta = 30, phi = 30, expand = 0.5, col = 8, main=title)
    }
    
    #write all the point info to a file
    info_file = sub("archv.tif", "info_full.csv",fixfile)
    write.csv(info, file=info_file, row.names = F) 
    
    dev.off() #turn off the plotting device
    
    #subset the points that passed the surface tests
    these = which(info[,"decision"] == 1)
    
    #if the number of sample points is less than 10 delete the image and return
    if(length(these) < 10){
      delete_files(fixfile, 2)
      return(0)
    } else {
      
      info = info[these,]
      
      #filter points based on rmse contribution
      maxr = endit = 10
      while(maxr >2 & endit != 0){
        xresid = (info[,"refx"]-info[,"fixx"])^2 #get the residuals of each x
        yresid = (info[,"refy"]-info[,"fixy"])^2 #get the residuals of each y
        r = (sqrt(xresid+yresid))/reso #get the rmse of each xy point
        totx = (1/length(info[,"refx"]))*(sum(xresid)) #intermediate step 
        toty = (1/length(info[,"refy"]))*(sum(yresid)) #intermediate step
        tot = sqrt(totx+toty)/reso #total rmse including all points
        if (tot != 0){contr = r/tot} else contr = r #error contribution of each point
        maxr = max(contr) #while loop controler
        b = which(contr < 2) #subset finder - is point 2 times or greater in contribution
        info = info[b,] #subset the info based on good rsme
        endit = sum(contr[b]) #while loop controler
      }
      
      #write out the filtered points that will be used in the transformation
      info_file = sub("archv.tif", "info_sub.csv",fixfile)
      write.csv(info, file=info_file, row.names = F)
      
      #adjust so that the coord is center of pixel
      info[,"fixx"] = info[,"fixx"]+(reso/2)
      info[,"fixy"] = info[,"fixy"]-(reso/2)
      
      #get the projection from the fix image
      wktfile = sub("archv.tif","wkt.txt", fixfile)
      projcmd = paste("gdalsrsinfo -o wkt", fixfile)
      proj = system(projcmd, intern = TRUE)
      write(proj, wktfile)
      
      #get the gcp string made
      fixcol = paste(info[,"fixcol"]) #fix col for tie point
      fixrow = paste(info[,"fixrow"]) #fix row for tie point
      refx = paste(info[,"refx"])  #fix x tie point coord 
      refy = paste(info[,"refy"])  #fix y tie point coord
      gcpstr = paste(" -gcp", fixcol, fixrow, refx, refy, collapse="")
      gcpfile = sub("archv.tif", "gcp.txt", fixfile)
      write(paste("reference file =", reffile), file=gcpfile)
      write(gcpstr, file=gcpfile, append=T)
      
      #gdal translate command
      tempname = sub("archv", "temp", fixfile) #"K:/scenes/034032/images/1976/LM10360321976248_archv.tif" 
      gdaltrans_cmd = paste("gdal_translate -of Gtiff -ot Byte -co INTERLEAVE=BAND -a_srs", wktfile, fixfile, tempname, gcpstr)
      system(gdaltrans_cmd)
                  
      #gdal warp command
      gdalwarp_cmd = paste("gdalwarp -of Gtiff -order 2 -ot Byte -srcnodata 0 -dstnodata 0 -co INTERLEAVE=BAND -overwrite -multi -tr", reso, reso, tempname, fixfile) #fixfile   "-tps"  "-order 2", "-order 3" 
      system(gdalwarp_cmd)
      
      #delete the temp file
      unlink(list.files(dirname(fixfile), pattern = "temp", full.names = T))
      return(1)
    } 
  } else {return(0)}
}
  
  
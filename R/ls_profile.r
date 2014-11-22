#' plot temporal spectral profile for points
#'
#' plot temporal spectral profile for points
#' @param msswrs1dir character. mss wrs1 directory path
#' @param msswrs2dir character. mss wrs2 directory path
#' @param tmwrs2dir character. tm wrs2 directory path
#' @param index character. spectral index to plot temporal profile
#' @param coords matrix or data.frame. matrix or dataframe with two row and n columns, column one id x coord, column two is y coord (optional)
#' @param n_random numeric. if coords are not given how many random points to plot temporal spectral profile
#' @param output character. how to display the profile(s). options "viewer" or "pdf"
#' @param mode. character. what to plot. options "apply" or "evaluate"
#' @import raster
#' @import ggplot2
#' @import zoo
#' @import segmented
#' @export

ls_profile = function(msswrs1dir, msswrs2dir, tmwrs2dir, index="tca", coords=NULL, n_random=10, output="viewer", mode="apply"){
  if(index == "tca"){search="tca.tif"}
  if(index == "tcb"){search="tc.tif"}
  if(index == "tcg"){search="tc.tif"}
  if(index == "tcw"){search="tc.tif"}
  
  msswrs1files = list.files(msswrs1dir, search, recursive=T, full.names=T)
  msswrs2files = list.files(msswrs2dir, search, recursive=T, full.names=T)
  tmwrs2files = list.files(tmwrs2dir, search, recursive=T, full.names=T)
  
  imgfiles = c(msswrs1files,msswrs2files,tmwrs2files)
  cloudfiles = sub(search, "cloudmask.tif",imgfiles)
  
  int = get_intersection(cloudfiles)
  img = crop(raster(imgfiles[1]),int)
  if(is.null(coords)){
    coords = sampleRandom(img, size=n_random, na.rm=TRUE, xy=T)
  }
  coords = data.frame(point=seq(1:nrow(coords)),x=coords[,1],y=coords[,2])
  
  if(index == "tca"){band=1}
  if(index == "tcb"){band=1}
  if(index == "tcg"){band=2}
  if(index == "tcw"){band=3}
  
  for(f in 1:length(imgfiles)){  #
    print(paste(f,"/",length(imgfiles),sep=""))
    bname = basename(imgfiles[f]) 
    img = raster(imgfiles[f],band=band)
    mask = raster(cloudfiles[f])
    value = extract(img, coords[,2:3])
    mvalue = extract(mask, coords[,2:3])
    these = which(mvalue == 0 | mvalue == NA)
    value[these] = NA
    yd = substr(bname,10,16)
    sat = substr(bname,1,3)
    if(substr(sat,2,2)=="M"){sensor = "MSS"}
    if(substr(sat,2,2)=="T"){sensor = "TM"}
    if(substr(sat,2,2)=="E"){sensor = "ETM+"}
    if(f == 1){fulldf = data.frame(coords,value,yd,sat,sensor)} else {
      df = data.frame(coords,value,yd,sat,sensor)
      fulldf = rbind(fulldf, df)
    }
  }
  
  fulldf$year = as.numeric(substr(fulldf$yd,1,4))
  
  if(output=="pdf"){pdf(outfile, width=11, height=5.5)}
  
  for(i in 1:length(unique(fulldf$point))){
    
    dfsub = subset(fulldf, point == i)
    
    #find the offset
    dfsubmss = subset(dfsub, sensor == "MSS")
    dfsubtm = subset(dfsub, sensor == "TM")
    
    thesemss = which(dfsubmss$yd %in% dfsubtm$yd)
    mss = dfsubmss[thesemss,]
    
    thesetm = which(dfsubtm$yd %in% dfsubmss$yd)
    tm = dfsubtm[thesetm,]
    
    if(length(mss) == length(tm)){
      mss = mss[order(mss$yd),]
      tm = tm[order(tm$yd),]
    } else {stop}
    
    offset = round(mean(mss$value - tm$value, na.rm=T)) #subtract this from all MSS
    
    #adjust the mss values
    dfsubmss$value = dfsubmss$value - offset
    dfsubmss$sensor = "MSS adj"
    dfsubadj = rbind(dfsub, dfsubmss)
    
    med = aggregate(value ~ year, data = dfsub, median, na.rm=T)
    mylm = lm(value ~ year, data = med)
    myseg = 6
    while (class(myseg) == "numeric"){
      myseg = ls_segit(mylm, seg.Z = ~ year, psi = NA, seg.control(stop.if.error = F, K = myseg, n.boot=0,it.max=20))
    }
    myfitted <- round(fitted(myseg))
    med$fitted = myfitted
    
    medadj = aggregate(value ~ year, data = dfsubadj, median, na.rm=T)
    mylm = lm(value ~ year, data = medadj)
    myseg = 6
    while (class(myseg) == "numeric"){
      myseg = ls_segit(mylm, seg.Z = ~ year, psi = NA, seg.control(stop.if.error = F, K = myseg, n.boot=0,it.max=20))
    }
    myfitted <- round(fitted(myseg))
    medadj$fitted = myfitted
    
    if(index == "tcb"){limits = c(-500, 10000)}
    if(index == "tcg"){limits = c(-500, 5000)}
    if(index == "tcw"){limits = c(-6000, 1000)}
    if(index == "tca"){limits = c(-500, 5000)}
    
    if(mode == "evaluate"){
      g= ggplot() +
        #geom_point(data=dfsubadj, aes(x=year, y=value, color=sensor), size = 4) + #, position = "jitter"
        ggtitle(paste("x =",dfsubadj$x[1],"y =",dfsubadj$y[1], "MSS mean offset =", offset))+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 90)) +
        geom_line(data=medadj, aes(x=year, y=fitted, color="MSS adj"))+
        geom_point(data=medadj, aes(x=year, y=value, color="MSS adj"), size=4, alpha=.7) +
        geom_point(data=med, aes(x=year, y=value, color="MSS"), size=4, alpha=.7) +
        geom_line(data=med, aes(x=year, y=fitted, color="MSS")) +
        ylim(limits) +
        xlim(1972,2014)
      ylab(index)
    }
    if(mode == "apply"){
      g= ggplot() +
        #geom_point(data=dfsubadj, aes(x=year, y=value, color=sensor), size = 4) + #, position = "jitter"
        ggtitle(paste("x =",round(dfsubadj$x[1]),"y =",round(dfsubadj$y[1])))+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 90)) +
        geom_point(data=med, aes(x=year, y=value), size=4) +
        geom_line(data=med, aes(x=year, y=fitted)) +
        ylim(limits) +
        xlim(1972,2014)
      ylab(index)
    }
    print(g)
  }
  
  if(output=="pdf"){dev.off()}
}



#' Calibrate oli images to TM images
#'
#' Calibrate oli images to TM images using linear regression
#' @param oliwrs2dir character. oli WRS-2 scene directory path
#' @param tmwrs2dir character. TM WRS-2 scene directory path
#' @import raster
#' @import MASS
#' @export


olical_single = function(oli_file, tm_file, overwrite=F){
  
  get_intersection = function(files){
    int = intersect(extent(raster(files[1])),extent(raster(files[2])))
    if(length(files) >= 3){for(i in 3:length(files))int = intersect(extent(raster(files[i])), int)}
    return(int)
  }
  
  predict_oli_index = function(tbl, outsampfile){  
    
    #create a multivariable linear model
    model = rlm(refsamp ~ b2samp + b3samp + b4samp + b5samp + b6samp + b7samp, data=tbl) #
    
    tbl$singlepred = round(predict(model))
    write.csv(tbl, outsampfile, row.names=F)
    
    #plot the regression
    r = cor(tbl$refsamp, tbl$singlepred)
    coef = rlm(tbl$refsamp ~ tbl$singlepred)
    
    pngout = sub("samp.csv", "plot.png",outsampfile)
    png(pngout,width=700, height=700)
    title = paste(tbl$index[1],"linear regression: slope =",paste(signif(coef$coefficients[2], digits=3),",",sep=""),
                  "y Intercept =",paste(round(coef$coefficients[1], digits=3),",",sep=""),
                  "r =",signif(r, digits=3))
    plot(x=tbl$singlepred,y=tbl$refsamp,
         main=title,
         xlab=paste(tbl$oli_img[1],tbl$index[1]),
         ylab=paste(tbl$ref_img[1],tbl$index[1]))
    abline(coef = coef$coefficients, col="red")  
    dev.off()
    
    #return the information
    coef_tbl = data.frame(rbind(model$coefficients))
    cnames = c("yint","b2c","b3c","b4c","b5c","b6c","b7c")
    colnames(coef_tbl) = cnames
    tbls = list(coef_tbl,tbl)
    return(tbls)
  }

  #define the filenames
  oli_sr_file = oli_file
  oli_mask_file = sub("l8sr.tif", "cloudmask.tif", oli_sr_file)
  ref_tc_file = tm_file
  ref_tca_file = sub("tc", "tca", ref_tc_file)
  ref_mask_file = sub("tc", "cloudmask", ref_tc_file)
  
  #make new directory
  dname = dirname(oli_sr_file)
  oliimgid = substr(basename(oli_sr_file),1,16)
  outdir = file.path(substr(dname,1,nchar(dname)-12),"calibration", oliimgid)  #-5
  dir.create(outdir, showWarnings = F, recursive=T)
  
  #check to see if single cal has already been run
  files = list.files(outdir)
  thesefiles = c("tca_cal_plot.png","tcb_cal_plot.png","tcg_cal_plot.png","tcw_cal_plot.png",
                 "tca_cal_samp.csv","tcb_cal_samp.csv","tcg_cal_samp.csv","tcw_cal_samp.csv")
  results = rep(NA,length(thesefiles))
  for(i in 1:length(results)){
    test = grep(thesefiles[i], files)
    results[i] = length(test) > 0
  }
  if(all(results) == T & overwrite == F){return(0)}
  
  
  #load files as raster
  oli_sr_img = brick(oli_sr_file)
  oli_mask_img = raster(oli_mask_file)
  ref_tc_img = brick(ref_tc_file)
  ref_tca_img  = raster(ref_tca_file)
  ref_mask_img = raster(ref_mask_file)
  
  #align the extents
  extent(oli_sr_img)  = alignExtent(oli_sr_img, ref_tc_img, snap="near")
  extent(oli_mask_img) = alignExtent(oli_mask_img, ref_tc_img, snap="near")
  extent(ref_tc_img)   = alignExtent(ref_tc_img, ref_tc_img, snap="near")
  extent(ref_tca_img)  = alignExtent(ref_tca_img, ref_tc_img, snap="near")
  extent(ref_mask_img) = alignExtent(ref_mask_img, ref_tc_img, snap="near")
  
  #crop the images to their intersection
  int = get_intersection(c(oli_mask_file,ref_mask_file))
  oli_b5_img = crop(subset(oli_sr_img,5),int)
  ref_tca_img = crop(ref_tca_img,int)
  oli_mask_img = crop(oli_mask_img,int)
  ref_mask_img = crop(ref_mask_img,int)
  
  #make a composite mask

  oli_mask_v = as.vector(oli_mask_img)
  ref_mask_v = as.vector(ref_mask_img)

  mask = oli_mask_v*ref_mask_v #make composite mask
  oli_mask_v = ref_mask_v = 0 # save memory
  
  #load oli and etm+ bands
  oli_b5_v = as.vector(oli_b5_img)
  ref_tca_v = as.vector(ref_tca_img)
  
  dif = oli_b5_v - ref_tca_v #find the difference
  oli_b5_v = ref_tca_v = 0 #save memory
  nas = which(mask == 0) #find the bads in the mask
  dif[nas] = NA #set the bads in the dif to NA so they are not included in the calc of mean and stdev
  stdv = sd(dif, na.rm=T) #get stdev of difference
  center = mean(dif, na.rm=T) #get the mean difference
  dif = dif < (center+stdv*2) & dif > (center-stdv*2) #find the pixels that are not that different
    
  
  goods = which(dif == 1)
  if(length(goods) < 20000){return(0)}
  
  #random sample
  samp = sample(1:length(goods), 20000)
  samp = goods[samp]
  sampxy = xyFromCell(oli_mask_img, samp)
  
  #save memory
  mask = 0
  
  #extract the sample pixels from the bands
  olisamp = extract(subset(oli_sr_img, 2:7), sampxy)
  tcsamp = extract(ref_tc_img, sampxy)
  tcasamp = extract(ref_tca_img, sampxy)
  
  #make sure the values are good for running regression on (diversity)
  unib2samp = length(unique(olisamp[,1]))
  unib3samp = length(unique(olisamp[,2]))
  unib4samp = length(unique(olisamp[,3]))
  unib5samp = length(unique(olisamp[,4]))
  unib6samp = length(unique(olisamp[,5]))
  unib7samp = length(unique(olisamp[,6]))
  
  unitcbsamp = length(unique(tcsamp[,1]))
  unitcgsamp = length(unique(tcsamp[,2]))
  unitcwsamp = length(unique(tcsamp[,3]))
  unitcasamp = length(unique(tcasamp))
  
  
  if(unib2samp < 15 | unib3samp < 15 | unib4samp < 15 | unib5samp < 15 | unib6samp < 15 | 
     unib7samp < 15 | unitcbsamp < 15 | unitcgsamp < 15 | unitcwsamp < 15 | unitcasamp < 15){return()}
  
  olibname = basename(oli_sr_file)
  refbname = basename(ref_tc_file)
  refabname = basename(ref_tca_file)
  
  tcb_tbl = data.frame(olibname,refbname,"tcb",sampxy,tcsamp[,1],olisamp)
  tcg_tbl = data.frame(olibname,refbname,"tcg",sampxy,tcsamp[,2],olisamp)
  tcw_tbl = data.frame(olibname,refbname,"tcw",sampxy,tcsamp[,3],olisamp)
  tca_tbl = data.frame(olibname,refabname,"tca",sampxy,tcasamp,olisamp)
  
  tcb_tbl = tcb_tbl[complete.cases(tcb_tbl),]
  tcg_tbl = tcg_tbl[complete.cases(tcg_tbl),]
  tcw_tbl = tcw_tbl[complete.cases(tcw_tbl),]
  tca_tbl = tca_tbl[complete.cases(tca_tbl),]
  
  cnames = c("oli_img","ref_img","index","x","y","refsamp","b2samp","b3samp","b4samp","b5samp","b6samp","b7samp") 
  colnames(tcb_tbl) = cnames
  colnames(tcg_tbl) = cnames
  colnames(tcw_tbl) = cnames
  colnames(tca_tbl) = cnames
  
  #predict the indices
  #TCB
  outsampfile = file.path(outdir,paste(oliimgid,"_tcb_cal_samp.csv",sep=""))
  model = predict_oli_index(tcb_tbl, outsampfile)
  bcoef = model[[1]]
  bsamp = model[[2]]
  br = cor(bsamp$refsamp, bsamp$singlepred)
  
  #TCG
  outsampfile = file.path(outdir,paste(oliimgid,"_tcg_cal_samp.csv",sep=""))
  model = predict_oli_index(tcg_tbl, outsampfile)
  gcoef = model[[1]]
  gsamp = model[[2]]
  gr = cor(gsamp$refsamp, gsamp$singlepred)
  
  #TCW
  outsampfile = file.path(outdir,paste(oliimgid,"_tcw_cal_samp.csv",sep=""))
  model = predict_oli_index(tcw_tbl, outsampfile)
  wcoef = model[[1]]
  wsamp = model[[2]]
  wr = cor(wsamp$refsamp, wsamp$singlepred)
  
  #TCA
  outsampfile = file.path(outdir,paste(oliimgid,"_tca_cal_samp.csv",sep=""))
  model = predict_oli_index(tca_tbl, outsampfile)
  acoef = model[[1]]
  asamp = model[[2]]
  ar = cor(asamp$refsamp, asamp$singlepred)
  
  #write out the coef files
  tcbinfo = data.frame(oli_file=olibname, ref_file=refbname, index="tcb", bcoef, r=br)
  tcginfo = data.frame(oli_file=olibname, ref_file=refbname, index="tcg", gcoef, r=gr)
  tcwinfo = data.frame(oli_file=olibname, ref_file=refbname, index="tcw", wcoef, r=wr)
  tcainfo = data.frame(oli_file=olibname, ref_file=refabname, index="tca", acoef, r=ar)
  
  tcbcoefoutfile = file.path(outdir,paste(oliimgid,"_tcb_cal_coef.csv",sep=""))
  tcgcoefoutfile = file.path(outdir,paste(oliimgid,"_tcg_cal_coef.csv",sep=""))
  tcwcoefoutfile = file.path(outdir,paste(oliimgid,"_tcw_cal_coef.csv",sep=""))
  tcacoefoutfile = file.path(outdir,paste(oliimgid,"_tca_cal_coef.csv",sep=""))
  
  write.csv(tcbinfo, tcbcoefoutfile, row.names=F)
  write.csv(tcginfo, tcgcoefoutfile, row.names=F)
  write.csv(tcwinfo, tcwcoefoutfile, row.names=F)
  write.csv(tcainfo, tcacoefoutfile, row.names=F)
}

#' Calibrate MSS images to TM images
#'
#' Calibrate MSS images to TM images using linear regression
#' @param msswrs2dir character. MSS WRS-2 scene directory path
#' @param tmwrs2dir character. TM WRS-2 scene directory path
#' @import raster
#' @import MASS
#' @export


msscal_single = function(mss_file, tm_file){
  
  get_intersection = function(files){
    int = intersect(extent(raster(files[1])),extent(raster(files[2])))
    if(length(files) >= 3){for(i in 3:length(files))int = intersect(extent(raster(files[i])), int)}
    return(int)
  }
  
  predict_mss_index = function(tbl, outsampfile){  
    #create a multivariable linear model
    model = rlm(refsamp ~ b1samp + b2samp + b3samp + b4samp, data=tbl) #tbl replaced final 1/22/2016
    
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
    plot(x=tbl$singlepred,y=tbl$refsamp, #tbl replaced final 1/22/2016
         main=title,
         xlab=paste(tbl$mss_img[1],tbl$index[1]),
         ylab=paste(tbl$ref_img[1],tbl$index[1]))   
    abline(coef = coef$coefficients, col="red")  
    dev.off()
    
    coef_tbl = data.frame(rbind(model$coefficients))
    cnames = c("yint","b1c","b2c","b3c","b4c")
    colnames(coef_tbl) = cnames
    tbls = list(coef_tbl,tbl)
    return(tbls)
  }
  
  #define the filenames
  mss_sr_file = mss_file
  mss_mask_file = sub("dos_sr_30m.tif", "cloudmask_30m.tif", mss_sr_file)
  ref_tc_file = tm_file
  ref_tca_file = sub("tc", "tca", ref_tc_file)
  ref_mask_file = sub("tc", "cloudmask", ref_tc_file)
  
  #make new directory
  dname = dirname(mss_sr_file)
  mssimgid = substr(basename(mss_sr_file),1,16)
  outdir = file.path(substr(dname,1,nchar(dname)-12),"calibration", mssimgid)  #-5
  dir.create(outdir, showWarnings = F, recursive=T)
  
  #load files as raster
  mss_sr_img = brick(mss_sr_file)
  mss_mask_img = raster(mss_mask_file)
  ref_tc_img = brick(ref_tc_file)
  ref_tca_img  = raster(ref_tca_file)
  ref_mask_img = raster(ref_mask_file)
  
  #align the extents
  extent(mss_sr_img)  = alignExtent(mss_sr_img, ref_tc_img, snap="near")
  extent(mss_mask_img) = alignExtent(mss_mask_img, ref_tc_img, snap="near")
  extent(ref_tc_img)   = alignExtent(ref_tc_img, ref_tc_img, snap="near")
  extent(ref_tca_img)  = alignExtent(ref_tca_img, ref_tc_img, snap="near")
  extent(ref_mask_img) = alignExtent(ref_mask_img, ref_tc_img, snap="near")
  
  #crop the images to their intersection
  int = get_intersection(c(mss_sr_file,mss_mask_file,ref_tc_file,ref_tca_file,ref_mask_file))
  mss_sr_img = crop(mss_sr_img,int)
  mss_mask_img = crop(mss_mask_img,int)
  ref_tc_img = crop(ref_tc_img,int)
  ref_tca_img = crop(ref_tca_img,int)
  ref_mask_img = crop(ref_mask_img,int)
  
  #make a composite mask
  mss_mask_v = as.vector(mss_mask_img)
  ref_mask_v = as.vector(ref_mask_img)
  mask = mss_mask_v*ref_mask_v
  mss_mask_v = ref_mask_v = 0 # save memory
  
  goods = which(mask == 1)
  if(length(goods) < 20000){return()}
  
  #random sample
  samp = sample(1:length(goods), 20000)
  samp = goods[samp]
  sampxy = xyFromCell(mss_mask_img, samp) #added on 1/22/2016
  
  #save memory
  mask = 0
  
  msssamp = extract(mss_sr_img, sampxy)
  tcsamp = extract(ref_tc_img, sampxy)
  tcasamp = extract(ref_tca_img, sampxy)
  
  #make sure the values are good for running regression on (diversity)
  unib1samp = length(unique(msssamp[,1]))
  unib2samp = length(unique(msssamp[,2]))
  unib3samp = length(unique(msssamp[,3]))
  unib4samp = length(unique(msssamp[,4]))
  
  unitcbsamp = length(unique(tcsamp[,1]))
  unitcgsamp = length(unique(tcsamp[,2]))
  unitcwsamp = length(unique(tcsamp[,3]))
  unitcasamp = length(unique(tcasamp))
  
  
  #if(unib1samp < 15 | unib2samp < 15 | unib3samp < 15 | unib4samp < 15 ){return()}
  if(unib1samp < 15 | unib2samp < 15 | unib3samp < 15 | unib4samp < 15 |
     unitcbsamp < 15 | unitcgsamp < 15 | unitcwsamp < 15 | unitcasamp < 15){return()}
  
  mssbname = basename(mss_sr_file)
  refbname = basename(ref_tc_file)
  refabname = basename(ref_tca_file)
  
  tcb_tbl = data.frame(mssbname,refbname,"tcb",sampxy,tcsamp[,1],msssamp)
  tcg_tbl = data.frame(mssbname,refbname,"tcg",sampxy,tcsamp[,2],msssamp)
  tcw_tbl = data.frame(mssbname,refbname,"tcw",sampxy,tcsamp[,3],msssamp)
  tca_tbl = data.frame(mssbname,refabname,"tca",sampxy,tcasamp,msssamp)
  
  tcb_tbl = tcb_tbl[complete.cases(tcb_tbl),]
  tcg_tbl = tcg_tbl[complete.cases(tcg_tbl),]
  tcw_tbl = tcw_tbl[complete.cases(tcw_tbl),]
  tca_tbl = tca_tbl[complete.cases(tca_tbl),]
  
  cnames = c("mss_img","ref_img","index","x","y","refsamp","b1samp","b2samp","b3samp","b4samp") 
  colnames(tcb_tbl) = cnames
  colnames(tcg_tbl) = cnames
  colnames(tcw_tbl) = cnames
  colnames(tca_tbl) = cnames
  
  outsampfile = file.path(outdir,paste(mssimgid,"_tcb_cal_samp.csv",sep=""))
  model = predict_mss_index(tcb_tbl, outsampfile)
  bcoef = model[[1]]
  bsamp = model[[2]]
  br = cor(bsamp$refsamp, bsamp$singlepred)
  
  outsampfile = file.path(outdir,paste(mssimgid,"_tcg_cal_samp.csv",sep=""))
  model = predict_mss_index(tcg_tbl, outsampfile)
  gcoef = model[[1]]
  gsamp = model[[2]]
  gr = cor(gsamp$refsamp, gsamp$singlepred)
  
  outsampfile = file.path(outdir,paste(mssimgid,"_tcw_cal_samp.csv",sep=""))
  model = predict_mss_index(tcw_tbl, outsampfile)
  wcoef = model[[1]]
  wsamp = model[[2]]
  wr = cor(wsamp$refsamp, wsamp$singlepred)
  
  #TCA
  outsampfile = file.path(outdir,paste(mssimgid,"_tca_cal_samp.csv",sep=""))
  model = predict_mss_index(tca_tbl, outsampfile)
  acoef = model[[1]]
  asamp = model[[2]]
  ar = cor(asamp$refsamp, asamp$singlepred)
  
  tcbinfo = data.frame(mss_file=mssbname, ref_file=refbname, index="tcb", bcoef, r=br)
  tcginfo = data.frame(mss_file=mssbname, ref_file=refbname, index="tcg", gcoef, r=gr)
  tcwinfo = data.frame(mss_file=mssbname, ref_file=refbname, index="tcw", wcoef, r=wr)
  tcainfo = data.frame(mss_file=mssbname, ref_file=refabname, index="tca", acoef, r=ar)
  
  tcbcoefoutfile = file.path(outdir,paste(mssimgid,"_tcb_cal_coef.csv",sep=""))
  tcgcoefoutfile = file.path(outdir,paste(mssimgid,"_tcg_cal_coef.csv",sep=""))
  tcwcoefoutfile = file.path(outdir,paste(mssimgid,"_tcw_cal_coef.csv",sep=""))
  tcacoefoutfile = file.path(outdir,paste(mssimgid,"_tca_cal_coef.csv",sep=""))
  
  write.csv(tcbinfo, tcbcoefoutfile, row.names=F)
  write.csv(tcginfo, tcgcoefoutfile, row.names=F)
  write.csv(tcwinfo, tcwcoefoutfile, row.names=F)
  write.csv(tcainfo, tcacoefoutfile, row.names=F)
}

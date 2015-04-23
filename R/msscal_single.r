#' Calibrate MSS images to TM images
#'
#' Calibrate MSS images to TM images using linear regression
#' @param msswrs2dir character. MSS WRS-2 scene directory path
#' @param tmwrs2dir character. TM WRS-2 scene directory path
#' @import raster
#' @import ggplot2
#' @import gridExtra
#' @export


msscal_single = function(mss_file, tm_file){
  
#   mss_file = mssf[1]
#   tm_file = tmf[1]
  
  get_intersection = function(files){
    int = intersect(extent(raster(files[1])),extent(raster(files[2])))
    if(length(files) >= 3){for(i in 3:length(files))int = intersect(extent(raster(files[i])), int)}
    return(int)
  }
  
  write_coef = function(mss_file, ref_file, index, coef,r){
    info = data.frame(mss_file = basename(mss_file),
                      ref_file = basename(ref_file),
                      index = index,
                      yint = as.numeric(coef[1]),
                      b1c = as.numeric(coef[2]),
                      b2c = as.numeric(coef[3]),
                      b3c = as.numeric(coef[4]),
                      b4c = as.numeric(coef[5]),
                      r=r)
    
    coefoutfile = file.path(outdir,paste(mssimgid,"_",index,"_cal_coef.csv",sep=""))
    write.csv(info, coefoutfile, row.names=F)
  }
  
  sample_it = function(img, bins, n){
    
    mi = min(img, na.rm=T)
    ma = max(img, na.rm=T)
    
    step = (ma - mi)/bins
    breaks = seq(mi,ma,step)
    
    min_samp = array(n, bins)
    for(i in 1:(length(breaks)-1)){
      these = which(img > breaks[i] & img <= breaks[i+1])
      if(i == 1){samp = sample(these, size=min(min_samp[i],length(these)))} else {
        samp = c(samp, sample(these, size=min(min_samp[i],length(these))))
      } 
    }
    return(samp)
  }
  
  #define the filenames
  mss_sr_file = mss_file
  mss_mask_file = sub("dos_sr_30m.tif", "cloudmask_30m.tif", mss_sr_file)
  ref_tc_file = tm_file
  ref_tca_file = sub("tc", "tca", ref_tc_file)
  ref_mask_file = sub("tc", "cloudmask", ref_tc_file)
  
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
  mss_mask_img = as.matrix(mss_mask_img)
  ref_mask_img = as.matrix(ref_mask_img)
  mask = mss_mask_img*ref_mask_img
  goods = which(mask == 1)
  refpix = as.matrix(ref_tca_img)[goods]
  
  #samp = sample_it(refpix, bins=20, n=1000)
  samp = sample(1:length(goods), 20000)
  samp = goods[samp]
  
  mss_mask_img = ref_mask_img = mask = refpix =  0
  
  b1samp = as.matrix(subset(mss_sr_img, 1))[samp]
  b2samp = as.matrix(subset(mss_sr_img, 2))[samp]
  b3samp = as.matrix(subset(mss_sr_img, 3))[samp]
  b4samp = as.matrix(subset(mss_sr_img, 4))[samp]
  samplen = length(samp)
  
  
  #predict the indices
  dname = dirname(mss_sr_file)
  mssimgid = substr(basename(mss_sr_file),1,16)
  outdir = file.path(substr(dname,1,nchar(dname)-12),"calibration", mssimgid)  #-5
  dir.create(outdir, showWarnings = F, recursive=T)
  
  #TCB
  refsamp = as.matrix(subset(ref_tc_img, 1))[samp]
  sampoutfile = file.path(outdir,paste(mssimgid,"_tcb_cal_samp.csv",sep=""))
  model = predict_mss_index(refsamp, b1samp, b2samp, b3samp, b4samp, mss_sr_file, ref_tc_file, "tcb", sampoutfile, samplen)
  bcoef = model[[1]]
  bsamp = model[[2]]
  #bplot = model[[3]]
  br = cor(bsamp$refsamp, bsamp$singlepred)
  
  #TCG
  refsamp = as.matrix(subset(ref_tc_img, 2))[samp]
  sampoutfile = file.path(outdir,paste(mssimgid,"_tcg_cal_samp.csv",sep=""))
  model = predict_mss_index(refsamp, b1samp, b2samp, b3samp, b4samp, mss_sr_file, ref_tc_file, "tcg", sampoutfile, samplen)
  gcoef = model[[1]]
  gsamp = model[[2]]
  #gplot = model[[3]]
  gr = cor(gsamp$refsamp, gsamp$singlepred)
  
  #TCW
  refsamp = as.matrix(subset(ref_tc_img, 3))[samp]
  sampoutfile = file.path(outdir,paste(mssimgid,"_tcw_cal_samp.csv",sep=""))
  model = predict_mss_index(refsamp, b1samp, b2samp, b3samp, b4samp, mss_sr_file, ref_tc_file, "tcw", sampoutfile, samplen)
  wcoef = model[[1]]
  wsamp = model[[2]]
  #wplot = model[[3]]
  wr = cor(wsamp$refsamp, wsamp$singlepred)
  
  #TCA
  singlepred = atan(gsamp$singlepred/bsamp$singlepred) * (180/pi) * 100
  refsamp = atan(gsamp$refsamp/bsamp$refsamp) * (180/pi) * 100
  tbl = data.frame(mss_img = rep(basename(mss_sr_file),length(singlepred)),
                   ref_img = rep(basename(ref_tc_file),length(singlepred)),
                   index = rep("tca",length(singlepred)),
                   refsamp,singlepred)
  final = tbl[complete.cases(tbl),]
  sampoutfile = file.path(outdir,paste(mssimgid,"_tca_cal_samp.csv",sep=""))
  write.csv(final, sampoutfile, row.names=F)
  
  
  r = cor(final$refsamp, final$singlepred)
  coef = rlm(final$refsamp ~ final$singlepred)
#   g = ggplot(final, aes(singlepred, refsamp)) +
#     stat_binhex(bins = 100)+
#     scale_fill_gradientn(name = "Count", colours = rainbow(7))+
#     xlab(paste(basename(mss_sr_file),"tca")) +
#     ylab(paste(basename(ref_tc_file),"tca")) +
#     ggtitle(paste("tca linear regression: slope =",paste(signif(coef$coefficients[2], digits=3),",",sep=""),
#                   "y Intercept =",paste(round(coef$coefficients[1], digits=3),",",sep=""),
#                   "r =",signif(r, digits=3))) +
#     theme(plot.title = element_text(size = 12)) +
#     geom_smooth(method="rlm", colour = "black", se = FALSE) + 
#     coord_fixed(ratio = 1)+
#     theme_bw()  
#   pngout = sub("samp.csv", "plot.png",sampoutfile)
#   png(pngout,width=700, height=700)
#   print(g)
#   dev.off()
  
  pngout = sub("samp.csv", "plot.png",sampoutfile)
  png(pngout,width=700, height=700)
  #print(g)
  title = paste("tca linear regression: slope =",paste(signif(coef$coefficients[2], digits=3),",",sep=""),
                                   "y Intercept =",paste(round(coef$coefficients[1], digits=3),",",sep=""),
                                   "r =",signif(r, digits=3))
  plot(x=final$singlepred,y=final$refsamp,
       main=title,
       xlab=paste(basename(mss_sr_file),"tca"),
       ylab=paste(basename(ref_tc_file),"tca"))
  abline(coef = coef$coefficients, col="red")  
  dev.off()
  
  info = data.frame(mss_file = basename(mss_sr_file), ref_file = basename(ref_tc_file),
                    index = "tca", yint = as.numeric(coef$coefficients[1]),
                    b1c = as.numeric(coef$coefficients[2]), r=r)
  
  coefoutfile = file.path(outdir,paste(mssimgid,"_tca_cal_coef.csv",sep=""))
  write.csv(info, coefoutfile, row.names=F)
  

  write_coef(mss_sr_file, ref_tc_file, "tcb", bcoef, br)
  write_coef(mss_sr_file, ref_tc_file, "tcg", gcoef, gr)
  write_coef(mss_sr_file, ref_tc_file, "tcw", wcoef, wr)

  
  #outfile = file.path(outdir,paste(mssimgid,"_tc_cal_planes.png",sep=""))
  #make_tc_planes_comparison(bsamp, gsamp, wsamp, outfile)
  
}

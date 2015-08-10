#' Calibrate oli images to TM images
#'
#' Calibrate oli images to TM images using linear regression
#' @param oliwrs2dir character. oli WRS-2 scene directory path
#' @param tmwrs2dir character. TM WRS-2 scene directory path
#' @import raster
#' @import ggplot2
#' @import gridExtra
#' @export


olical_single = function(oli_file, tm_file){
  
  #   oli_file = mssf[1]
  #   tm_file = tmf[1]
  
  get_intersection = function(files){
    int = intersect(extent(raster(files[1])),extent(raster(files[2])))
    if(length(files) >= 3){for(i in 3:length(files))int = intersect(extent(raster(files[i])), int)}
    return(int)
  }
  
  write_coef = function(oli_file, ref_file, index, coef,r){
    info = data.frame(oli_file = basename(oli_file),
                      ref_file = basename(ref_file),
                      index = index,
                      yint = as.numeric(coef[1]),
                      b1c = as.numeric(coef[2]),
                      b2c = as.numeric(coef[3]),
                      b3c = as.numeric(coef[4]),
                      b4c = as.numeric(coef[5]),
                      b5c = as.numeric(coef[6]),
                      b6c = as.numeric(coef[7]),
                      b7c = as.numeric(coef[8]),
                      r=r)
    
    coefoutfile = file.path(outdir,paste(mssimgid,"_",index,"_cal_coef.csv",sep=""))
    write.csv(info, coefoutfile, row.names=F)
  }
  
#   sample_it = function(img, bins, n){
#     
#     mi = min(img, na.rm=T)
#     ma = max(img, na.rm=T)
#     
#     step = (ma - mi)/bins
#     breaks = seq(mi,ma,step)
#     
#     min_samp = array(n, bins)
#     for(i in 1:(length(breaks)-1)){
#       these = which(img > breaks[i] & img <= breaks[i+1])
#       if(i == 1){samp = sample(these, size=min(min_samp[i],length(these)))} else {
#         samp = c(samp, sample(these, size=min(min_samp[i],length(these))))
#       } 
#     }
#     return(samp)
#   }
  
  #define the filenames
  oli_sr_file = oli_file
  oli_mask_file = sub("dos_sr_30m.tif", "cloudmask_30m.tif", oli_sr_file)
  ref_tc_file = tm_file
  ref_tca_file = sub("tc", "tca", ref_tc_file)
  ref_mask_file = sub("tc", "cloudmask", ref_tc_file)
  
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
  int = get_intersection(c(oli_sr_file,oli_mask_file,ref_tc_file,ref_tca_file,ref_mask_file))
  oli_sr_img = crop(oli_sr_img,int)
  oli_mask_img = crop(oli_mask_img,int)
  ref_tc_img = crop(ref_tc_img,int)
  ref_tca_img = crop(ref_tca_img,int)
  ref_mask_img = crop(ref_mask_img,int)
  
  #make a composite mask
  oli_mask_img = as.matrix(oli_mask_img)
  ref_mask_img = as.matrix(ref_mask_img)
  mask = oli_mask_img*ref_mask_img
  goods = which(mask == 1)
  if(length(goods) < 20000){return()}
  
  #stratified sample
  #refpix = as.matrix(ref_tca_img)[goods]
  #samp = sample_it(refpix, bins=20, n=1000)
  
  #random sample
  samp = sample(1:length(goods), 20000)
  samp = goods[samp]
  
  #save memory
  oli_mask_img = ref_mask_img = mask = refpix =  0
  
  #extract the sample pixels from the bands
  b1samp = as.matrix(subset(oli_sr_img, 1))[samp]
  b2samp = as.matrix(subset(oli_sr_img, 2))[samp]
  b3samp = as.matrix(subset(oli_sr_img, 3))[samp]
  b4samp = as.matrix(subset(oli_sr_img, 4))[samp]
  b5samp = as.matrix(subset(oli_sr_img, 5))[samp]
  b6samp = as.matrix(subset(oli_sr_img, 6))[samp]
  b7samp = as.matrix(subset(oli_sr_img, 7))[samp]
  
  #make sure the values are good for running regression on (diversity)
  unib1samp = length(unique(b1samp))
  unib2samp = length(unique(b2samp))
  unib3samp = length(unique(b3samp))
  unib4samp = length(unique(b4samp))
  unib5samp = length(unique(b5samp))
  unib6samp = length(unique(b6samp))
  unib7samp = length(unique(b7samp))
  if(unib1samp < 15 | unib2samp < 15 | unib3samp < 15 | unib4samp < 15 | unib5samp < 15 | unib6samp < 15 | unib7samp < 15){return()}
  
  samplen = length(samp)
  
  #predict the indices
  dname = dirname(oli_sr_file)
  oliimgid = substr(basename(oli_sr_file),1,16)
  outdir = file.path(substr(dname,1,nchar(dname)-12),"calibration", oliimgid)  #-5
  dir.create(outdir, showWarnings = F, recursive=T)
  
  #TCB
  refsamp = as.matrix(subset(ref_tc_img, 1))[samp]
  unirefsamp = length(unique(refsamp))
  if(unirefsamp < 15){return()}
  sampoutfile = file.path(outdir,paste(oliimgid,"_tcb_cal_samp.csv",sep=""))
  model = predict_oli_index(refsamp, b1samp, b2samp, b3samp, b4samp, b5samp, b6samp, b7samp, oli_sr_file, ref_tc_file, "tcb", sampoutfile, samplen)
  bcoef = model[[1]]
  bsamp = model[[2]]
  #bplot = model[[3]]
  br = cor(bsamp$refsamp, bsamp$singlepred)
  
  #TCG
  refsamp = as.matrix(subset(ref_tc_img, 2))[samp]
  unirefsamp = length(unique(refsamp))
  if(unirefsamp < 15){return()}
  sampoutfile = file.path(outdir,paste(oliimgid,"_tcg_cal_samp.csv",sep=""))
  model = predict_oli_index(refsamp, b1samp, b2samp, b3samp, b4samp, b5samp, b6samp, b7samp, oli_sr_file, ref_tc_file, "tcg", sampoutfile, samplen)
  gcoef = model[[1]]
  gsamp = model[[2]]
  #gplot = model[[3]]
  gr = cor(gsamp$refsamp, gsamp$singlepred)
  
  #TCW
  refsamp = as.matrix(subset(ref_tc_img, 3))[samp]
  unirefsamp = length(unique(refsamp))
  if(unirefsamp < 15){return()}
  sampoutfile = file.path(outdir,paste(oliimgid,"_tcw_cal_samp.csv",sep=""))
  model = predict_oli_index(refsamp, b1samp, b2samp, b3samp, b4samp, b5samp, b6samp, b7samp, oli_sr_file, ref_tc_file, "tcw", sampoutfile, samplen)
  wcoef = model[[1]]
  wsamp = model[[2]]
  #wplot = model[[3]]
  wr = cor(wsamp$refsamp, wsamp$singlepred)
  
  #TCA
  singlepred = atan(gsamp$singlepred/bsamp$singlepred) * (180/pi) * 100
  refsamp = atan(gsamp$refsamp/bsamp$refsamp) * (180/pi) * 100
  tbl = data.frame(oli_img = rep(basename(oli_sr_file),length(singlepred)),
                   ref_img = rep(basename(ref_tc_file),length(singlepred)),
                   index = rep("tca",length(singlepred)),
                   refsamp,singlepred)
  final = tbl[complete.cases(tbl),]
  sampoutfile = file.path(outdir,paste(oliimgid,"_tca_cal_samp.csv",sep=""))
  write.csv(final, sampoutfile, row.names=F)
  
  
  r = cor(final$refsamp, final$singlepred)
  coef = rlm(final$refsamp ~ final$singlepred)
  #   g = ggplot(final, aes(singlepred, refsamp)) +
  #     stat_binhex(bins = 100)+
  #     scale_fill_gradientn(name = "Count", colours = rainbow(7))+
  #     xlab(paste(basename(oli_sr_file),"tca")) +
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
       xlab=paste(basename(oli_sr_file),"tca"),
       ylab=paste(basename(ref_tc_file),"tca"))
  abline(coef = coef$coefficients, col="red")  
  dev.off()
  
  info = data.frame(oli_file = basename(oli_sr_file), ref_file = basename(ref_tc_file),
                    index = "tca", yint = as.numeric(coef$coefficients[1]),
                    b1c = as.numeric(coef$coefficients[2]), r=r)
  
  coefoutfile = file.path(outdir,paste(oliimgid,"_tca_cal_coef.csv",sep=""))
  write.csv(info, coefoutfile, row.names=F)
  
  
  write_coef(oli_sr_file, ref_tc_file, "tcb", bcoef, br)
  write_coef(oli_sr_file, ref_tc_file, "tcg", gcoef, gr)
  write_coef(oli_sr_file, ref_tc_file, "tcw", wcoef, wr)
  
  
  #outfile = file.path(outdir,paste(oliimgid,"_tc_cal_planes.png",sep=""))
  #make_tc_planes_comparison(bsamp, gsamp, wsamp, outfile)
  
}

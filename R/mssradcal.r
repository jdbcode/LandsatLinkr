#' Radiometrically calibrate MSS TC to TM TC
#'
#' Spectrally calibrate MSS TC to TM TC
#' @param mss_rad_file filename full path to MSS radiance file
#' @param ref_tc_file filename full path to TM TC reference file
#' @import raster
#' @export


mssradcal = function(mss_rad_file, ref_tc_file){
  bins=20
  n=10000
  
  mss_tc_file = sub("radiance.tif", "tc.tif", mss_rad_file)
  mss_tca_file = sub("radiance.tif", "tca.tif", mss_rad_file)
  mss_mask_file = sub("radiance.tif", "cloudmask.tif", mss_rad_file)
  ref_tca_file = sub("tc", "tca", ref_tc_file)
  ref_mask_file = sub("tc", "cloudmask", ref_tc_file)
  files = c(mss_rad_file, mss_tc_file, mss_tca_file, mss_mask_file,
            ref_tc_file, ref_tca_file, ref_mask_file)
  
  #load files as raster
  mss_rad_img = brick(files[1])
  mss_tc_img  = brick(files[2])
  mss_tca_img = raster(files[3])
  mss_mask_img = raster(files[4])
  ref_tc_img = brick(files[5])
  ref_tca_img  = raster(files[6])
  ref_mask_img = raster(files[7])
  
  #align the extents
  mss_rad_img_ex  = alignExtent(mss_rad_img, ref_tc_img, snap="near")
  mss_tc_img_ex   = alignExtent(mss_tc_img, ref_tc_img, snap="near")
  mss_tca_img_ex  = alignExtent(mss_tca_img, ref_tc_img, snap="near")
  mss_mask_img_ex = alignExtent(mss_mask_img, ref_tc_img, snap="near")
  ref_tc_img_ex   = alignExtent(ref_tc_img, ref_tc_img, snap="near")
  ref_tca_img_ex  = alignExtent(ref_tca_img, ref_tc_img, snap="near")
  ref_mask_img_ex = alignExtent(ref_mask_img, ref_tc_img, snap="near")
  
  extent(mss_rad_img) = mss_rad_img_ex
  extent(mss_tc_img) = mss_tc_img_ex
  extent(mss_tca_img) = mss_tca_img_ex
  extent(mss_mask_img) = mss_mask_img_ex
  extent(ref_tc_img) = ref_tc_img_ex
  extent(ref_tca_img) = ref_tca_img_ex
  extent(ref_mask_img) = ref_mask_img_ex
  
  #crop the images to their intersection
  int = get_intersection(files)
  mss_rad_img = crop(mss_rad_img,int)
  mss_tc_img = crop(mss_tc_img,int)
  mss_tca_img = crop(mss_tca_img,int)
  mss_mask_img = crop(mss_mask_img,int)
  ref_tc_img = crop(ref_tc_img,int)
  ref_tca_img = crop(ref_tca_img,int)
  ref_mask_img = crop(ref_mask_img,int)
  
  #make a composite mask
  mss_mask_img = as.matrix(mss_mask_img)
  ref_mask_img = as.matrix(ref_mask_img)
  mask = mss_mask_img*ref_mask_img #*fix_goods*ref_goods #make a composite mask
  zeros = which(mask == 0)
  mss_mask_img = ref_mask_img = mask = 0 #same memory
  
  #find points that have not changed drastically
  mss_tca_img = as.matrix(mss_tca_img)
  ref_tca_img = as.matrix(ref_tca_img)
  mss_tca_img[zeros] = NA
  ref_tca_img[zeros] = NA
  dif = mss_tca_img-ref_tca_img
  ave = mean(dif, na.rm=T)
  devi = sd(dif, na.rm=T)
  goods = which(dif < ave+devi & dif > ave-devi) #elements
  
  #create output names
  tcb = sub("tc.tif", "tcb.tif", mss_tc_file)
  tcg = sub("tc.tif", "tcg.tif", mss_tc_file)
  tcw = sub("tc.tif", "tcw.tif", mss_tc_file)
  tca= sub("tc.tif", "tca_cal.tif", mss_tc_file)
  tcall = sub("tc.tif", "tc_cal.tif", mss_tc_file)
  tcpdf = sub("tc.tif", "tc_cal.pdf", mss_tc_file)
  tcapdf = sub("tc.tif", "tca_cal.pdf", mss_tc_file)
  tccsv = sub("tc.tif", "tc_cal.csv", mss_tc_file)
  tcacsv = sub("tc.tif", "tca_cal.csv", mss_tc_file)
  
  outfile = c(tcb,tcg,tcw,tca,tcall,tcpdf,tcapdf,tccsv,tcacsv)
  
  ###############################################################
  #modeling
  #model tcb
  model = predict_mss_brt_grn_tca(subset(mss_tc_img, 1), subset(ref_tc_img, 1), goods, bins, n) 
  bcoef = model[[1]]
  bsamp = model[[2]]
  
  #model tcg
  model = predict_mss_brt_grn_tca(subset(mss_tc_img, 2), subset(ref_tc_img, 2), goods, bins, n) 
  gcoef = model[[1]]
  gsamp = model[[2]]
  
  #model tCw
  model = predict_mss_wetness(mss_rad_img, ref_tc_img, goods, bins, n)
  wcoef = model[[1]]
  wsamp = model[[2]]
  
  #model tca
  model = predict_mss_brt_grn_tca(mss_tca_img, ref_tca_img, goods, bins, n) 
  acoef = model[[1]]
  asamp = model[[2]]
  
  #write tcb
  write_tc(mss_tc_file, "tcb", bcoef, outfile[1])
  #write tcg
  write_tc(mss_tc_file, "tcg", gcoef, outfile[2])
  #write tcw
  write_tc(mss_rad_file, "tcw", wcoef, outfile[3])
  #extract info to make wetness scatterplot
  mss_tcw_img = raster(outfile[3])
  mss_tcw_img_ex = alignExtent(mss_tcw_img, ref_tc_img, snap="near")
  extent(mss_tcw_img) = mss_tcw_img_ex
  mss_tcw_img = crop(mss_tcw_img,int)
  mss_tcw_img = as.matrix(mss_tcw_img)
  mss_tcw_img = mss_tcw_img[goods]
  msssamp = mss_tcw_img[wsamp]
  ref_tcw_img = subset(ref_tc_img,3)
  ref_tcw_img = as.matrix(ref_tcw_img)
  ref_tcw_img = ref_tcw_img[goods]
  refsamp = ref_tcw_img[wsamp]
  wsamp = data.frame(refsamp,msssamp)
  model = rlm(refsamp ~ msssamp, data=wsamp)
  yint = model$coefficients[1]
  coef = model$coefficients[2]
  wcoef = data.frame(yint,coef)
  #write tca
  write_tc(mss_tca_file, "tca", acoef, outfile[4])
  
  #write out the tc stack
  tc_files = paste(outfile[1:3], collapse=" ")
  cmd = paste("gdal_merge.py -o",outfile[5],"-separate -a_nodata -32768", tc_files)
  shell(cmd)
  
  #calc correlation for plotting
  br = cor(bsamp$refsamp, bsamp$msssamp)
  gr = cor(gsamp$refsamp, gsamp$msssamp)
  wr = cor(wsamp$refsamp, wsamp$msssamp)
  ar = cor(asamp$refsamp, asamp$msssamp)
  
  #create scatterplots
  bplot = plot_it(bsamp, bcoef, br, outfile[5], ref_tc_file, "TC Brightness")
  gplot = plot_it(gsamp, gcoef, gr, outfile[5], ref_tc_file, "TC Greenness")
  wplot = plot_it(wsamp, wcoef, wr, outfile[5], ref_tc_file, "TC Wetness")
  aplot = plot_it(asamp, acoef, ar, outfile[4], ref_tca_file, "TC Angle")
  
  #save scatterplots as .pdf
  unlink(outfile[6])
  pdf(file = outfile[6], width=11, height=8.5)
  print(bplot)
  print(gplot)
  print(wplot)
  dev.off()
  
  unlink(outfile[7])
  pdf(file = outfile[7], width=11, height=8.5)
  print(aplot)
  dev.off()
  
  #make .csv summaries
  labels = c("index","ref_file","ref_mask_file","fix_file","fix_mask_file","slope","y_intercept","correlation")
  tcbinfo = data.frame("tc brightness",
                       basename(ref_tc_file),
                       basename(ref_mask_file),
                       basename(mss_tc_file),
                       basename(mss_mask_file),
                       bcoef$coef,
                       bcoef$yint,
                       br)
  colnames(tcbinfo) = labels
  tcginfo = data.frame("tc greenness",
                       basename(ref_tc_file),
                       basename(ref_mask_file),
                       basename(mss_tc_file),
                       basename(mss_mask_file),
                       gcoef$coef,
                       gcoef$yint,
                       gr)
  colnames(tcginfo) = labels
  tcwinfo = data.frame("tc wetness",
                       basename(ref_tc_file),
                       basename(ref_mask_file),
                       basename(mss_rad_file),
                       basename(mss_mask_file),
                       wcoef$coef,
                       wcoef$yint,
                       wr)
  colnames(tcwinfo) = labels
  tcainfo = data.frame("tc angle",
                       basename(ref_tca_file),
                       basename(ref_mask_file),
                       basename(mss_tca_file),
                       basename(mss_mask_file),
                       acoef$coef,
                       acoef$yint,
                       ar) 
  colnames(tcainfo) = labels
  tcbgw_info = rbind(tcbinfo,tcginfo,tcwinfo)
  write.csv(tcbgw_info, file=outfile[8], row.names = F)
  write.csv(tcainfo, file=outfile[9], row.names = F)
  unlink(outfile[1:3])
}

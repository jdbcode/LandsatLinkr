#' Plot the three TC planes for MSS and TM TC 
#'
#' Plot the three TC planes for MSS and TM TC. It should only be done for coincident images.
#' @param mss_tc_file character. full filename path to MSS tc file
#' @param tm_tc_file character. full filename path to TM TC file
#' @param bins integer. 
#' @param n integer.
#' @param outfile
#' @import raster
#' @import grid
#' @import gridExtra
#' @import ggplot2


plot_mss_tm_tc = function(mss_tc_file, tm_tc_file, bins, n, outfile){
  
  mss_tca_file = sub("tc_cal.tif", "tca_cal.tif", mss_tc_file)
  mss_mask_file = sub("tc_cal.tif", "cloudmask.tif", mss_tc_file)
  tm_tca_file = sub("tc.tif", "tca.tif.tif", tm_tc_file)
  tm_mask_file = sub("tc.tif", "cloudmask.tif", ref_tc_file)
  
  files = c(mss_tc_file, tm_tc_file, mss_mask_file, tm_mask_file, mss_tca_file, tm_tca_file)
  int = get_intersection(files)
  
  mss_tc_img = brick(mss_tc_file) 
  tm_tc_img = brick(tm_tc_file)
  mss_tca_img = raster(mss_tca_file)
  tm_tca_img = raster(tm_tca_file)
  mss_mask_img = raster(mss_mask_file)
  tm_mask_img = raster(tm_mask_file)
  
  mss_tc_imgex = alignExtent(mss_tc_img, tm_tc_img, snap="near") 
  mss_tca_imgex = alignExtent(mss_tca_img, tm_tc_img, snap="near")
  tm_tca_imgex = alignExtent(tm_tca_img, tm_tc_img, snap="near")
  mss_mask_imgex = alignExtent(mss_mask_img, tm_tc_img, snap="near") 
  tm_mask_imgex = alignExtent(tm_mask_img, tm_tc_img, snap="near")
  extent(mss_tc_img) = mss_tc_imgex
  extent(mss_tca_img) = mss_tca_imgex
  extent(tm_tca_img) = tm_tca_imgex
  extent(mss_mask_img) = mss_mask_imgex
  extent(tm_mask_img) = tm_mask_imgex
  mss_tc_img = crop(mss_tc_img,int)
  tm_tc_img = crop(tm_tc_img, int)
  mss_tca_img = crop(mss_tca_img,int)
  tm_tca_img = crop(tm_tca_img, int)
  mss_mask_img = crop(mss_mask_img,int)
  tm_mask_img = crop(tm_mask_img, int)
  
  #make a composite mask
  mss_mask_img = as.matrix(mss_mask_img)
  tm_mask_img = as.matrix(tm_mask_img)
  mask = mss_mask_img*tm_mask_img #*fix_goods*ref_goods #make a composite mask
  zeros = which(mask == 0)
  mss_mask_img = tm_mask_img = mask = 0 #same memory
  
  #find points that have not changed drastically
  mss_tca_img = as.matrix(mss_tca_img)
  tm_tca_img = as.matrix(tm_tca_img)
  mss_tca_img[zeros] = NA
  tm_tca_img[zeros] = NA
  dif = mss_tca_img-tm_tca_img
  ave = mean(dif, na.rm=T)
  devi = sd(dif, na.rm=T)
  goods = which(dif < ave+devi & dif > ave-devi) #elements
  mss_tca_img = tm_tca_img = 0
  
  #create samples for tc separately based on goods from above
  
  bsamp = sample_it(as.matrix(subset(tm_tc_img,1))[goods], bins, n)
  gsamp = sample_it(as.matrix(subset(tm_tc_img,2))[goods], bins, n)
  wsamp = sample_it(as.matrix(subset(tm_tc_img,3))[goods], bins, n)
  allsamp = unique(c(bsamp,gsamp,wsamp))
  
  #subset the samples from each tc index fpr each sensor 
  mss_tcb_img = as.matrix(subset(mss_tc_img,1))[goods][allsamp]
  mss_tcg_img = as.matrix(subset(mss_tc_img,2))[goods][allsamp]
  mss_tcw_img = as.matrix(subset(mss_tc_img,3))[goods][allsamp]
  tm_tcb_img = as.matrix(subset(tm_tc_img,1))[goods][allsamp]
  tm_tcg_img = as.matrix(subset(tm_tc_img,2))[goods][allsamp]
  tm_tcw_img = as.matrix(subset(tm_tc_img,3))[goods][allsamp]
  
  #make a big table
  df = data.frame(mss_tcb_img,
                  mss_tcg_img,
                  mss_tcw_img,
                  tm_tcb_img,
                  tm_tcg_img,
                  tm_tcw_img)
  
  blimits = c(-500, 10000)
  glimits = c(-500, 5000)
  wlimits = c(-6000, 1000)
  
  #create plots
  mssbg = ggplot(df, aes(mss_tcb_img, mss_tcg_img)) +
    stat_binhex(bins = 100)+
    scale_fill_gradientn(colours = rainbow(7)) +
    scale_x_continuous(limits = blimits) +
    scale_y_continuous(limits = glimits) +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  mssbw = ggplot(df, aes(mss_tcb_img, mss_tcw_img)) +
    stat_binhex(bins = 100)+
    scale_fill_gradientn(colours = rainbow(7)) +
    scale_x_continuous(limits = blimits) +
    scale_y_continuous(limits = wlimits) +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  msswg = ggplot(df, aes(mss_tcw_img, mss_tcg_img)) +
    stat_binhex(bins = 100)+
    scale_fill_gradientn(colours = rainbow(7)) +
    scale_x_continuous(limits = wlimits) +
    scale_y_continuous(limits = glimits) +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  tmbg = ggplot(df, aes(tm_tcb_img, tm_tcg_img)) +
    stat_binhex(bins = 100)+
    scale_fill_gradientn(colours = rainbow(7)) +
    scale_x_continuous(limits = blimits) +
    scale_y_continuous(limits = glimits) +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  tmbw = ggplot(df, aes(tm_tcb_img, tm_tcw_img)) +
    stat_binhex(bins = 100)+
    scale_fill_gradientn(colours = rainbow(7)) +
    scale_x_continuous(limits = blimits) +
    scale_y_continuous(limits = wlimits) +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  tmwg = ggplot(df, aes(tm_tcw_img, tm_tcg_img)) +
    stat_binhex(bins = 100)+
    scale_fill_gradientn(colours = rainbow(7)) +
    scale_x_continuous(limits = wlimits) +
    scale_y_continuous(limits = glimits) +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 90, hjust = 1))
  
  plot_list = c(list(mssbg),list(mssbw),list(msswg),list(tmbg),list(tmbw),list(tmwg))
  pdf(file = outfile, width=11, height=8.5)
  do.call(grid.arrange, c(plot_list, list(ncol=3)))
  dev.off() 
}
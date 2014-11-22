#' calibrate MSS TC brightness, greenness, and angle
#'
#' calibrate MSS TC brightness, greenness, and angle by linear regression between MSS TC components and corresponding reference TC components
#' @param mss_img raster* mss radiance brick or stack image containing bands 1-4 
#' @param ref_img raster* TC wetness reference image
#' @param goods array index of pixels to draw a sample from
#' @param bins integer number of histogram bins from which to sample
#' @param n integer number of pixels to sample per histogram bin


predict_mss_brt_grn_tca = function(mss_img, ref_img, goods, bins, n){
  #extract the goods
  mss_img = as.matrix(mss_img)[goods]
  ref_img = as.matrix(ref_img)[goods]
  
  #create a sample; 1 for brightness and one for greeness
  samp = sample_it(ref_img, bins, n)
  
  #extract the sample
  msssamp = mss_img[samp]
  refsamp = ref_img[samp]
  
  #make a dateframe of the training sample
  tbl = data.frame(refsamp, msssamp) #,tcbsamp,tcgsamp,tcysamp
  final = tbl[complete.cases(tbl),]
  
  #create a linear model
  model = rlm(refsamp ~ msssamp, data=final)
  yint = model$coefficients[1]
  coef = model$coefficients[2]
  coef = data.frame(yint,coef)  
  
  tbls = list(coef, final)
  return(tbls)
}
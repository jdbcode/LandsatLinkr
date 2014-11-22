#' predict mss TC wetness
#'
#' predicts tasseled cap wetness for mss by linear regression between mss bands 1-4 and TM wetness
#' @param mss_rad_img raster* mss radiance brick or stack image containing bands 1-4 
#' @param ref_tc_img raster* TC wetness reference image
#' @param goods array index of pixels to draw a sample from
#' @param bins integer number of histogram bins from which to sample
#' @param n integer number of pixels to sample per histogram bin
#' @import raster


predict_mss_wetness = function(mss_rad_img, ref_tc_img, goods, bins, n){
  #convert to matrix for faster processing
  b1 = as.matrix(subset(mss_rad_img, 1))[goods]
  b2 = as.matrix(subset(mss_rad_img, 2))[goods]
  b3 = as.matrix(subset(mss_rad_img, 3))[goods]
  b4 = as.matrix(subset(mss_rad_img, 4))[goods]
  ref_img = as.matrix(subset(ref_tc_img, 3))[goods]
  samp = sample_it(ref_img, bins, n)
  
  #extract the sample
  b1samp = b1[samp]
  b2samp = b2[samp]
  b3samp = b3[samp]
  b4samp = b4[samp]
  y = ref_img[samp]
  
  #make a dateframe of the training sample
  tbl = data.frame(y,b1samp,b2samp,b3samp,b4samp) #,tcbsamp,tcgsamp,tcysamp
  final = tbl[complete.cases(tbl),]
  
  #create a multivariable linear model
  model = rlm(y ~ b1samp + b2samp + b3samp + b4samp, data=final)
  
  #apply the model to the original img
  yint = model$coefficients[1]
  b1c = model$coefficients[2]
  b2c = model$coefficients[3]
  b3c = model$coefficients[4]
  b4c = model$coefficients[5]
  
  coef = data.frame(yint,b1c,b2c,b3c,b4c)
  tbls = list(coef,samp)
  return(tbls)
}
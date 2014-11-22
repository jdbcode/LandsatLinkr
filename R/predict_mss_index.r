#' Predict MSS TC index
#'
#' Predict MSS TC index
#' @param file The full path name of the *archv file
#' @param pngout Logical. If TRUE then a 4-panel .png image file will be created that for each MSS band displays a truncated histogram and the estimated Dark Object Value
#' @import MASS
#' @import ggplot2
#' @export


predict_mss_index = function(refsamp, b1samp, b2samp, b3samp, b4samp, mss_sr_file, ref_file, index, outsampfile, samplen){  
  #make a dateframe of the training sample
  tbl = data.frame(mss_img = rep(basename(mss_sr_file),samplen),
                   ref_img = rep(basename(ref_file),samplen),
                   index = rep(index,samplen),
                   refsamp,b1samp,b2samp,b3samp,b4samp)
  final = tbl[complete.cases(tbl),]
  
  
  #create a multivariable linear model
  model = rlm(refsamp ~ b1samp + b2samp + b3samp + b4samp, data=final)
  
  final$singlepred = round(predict(model))
  write.csv(final, outsampfile, row.names=F)
  
  #apply the model to the original img
  yint = model$coefficients[1]
  b1c = model$coefficients[2]
  b2c = model$coefficients[3]
  b3c = model$coefficients[4]
  b4c = model$coefficients[5]
  
  #plot the regression
  r = cor(final$refsamp, final$singlepred)
  coef = rlm(final$refsamp ~ final$singlepred)
  
  g = ggplot(final, aes(singlepred, refsamp)) +
    stat_binhex(bins = 100)+
    scale_fill_gradientn(name = "Count", colours = rainbow(7))+
    xlab(paste(basename(mss_sr_file),index)) +
    ylab(paste(basename(ref_file),index)) +
    ggtitle(paste(index,"linear regression: slope =",paste(signif(coef$coefficients[2], digits=3),",",sep=""),
                  "y Intercept =",paste(round(coef$coefficients[1], digits=3),",",sep=""),
                  "r =",signif(r, digits=3))) +
    theme(plot.title = element_text(size = 12)) +
    geom_smooth(method="rlm", colour = "black", se = FALSE) + 
    coord_fixed(ratio = 1)+
    theme_bw()
  
  pngout = sub("samp.csv", "plot.png",outsampfile)
  png(pngout,width=700, height=700)
  print(g)
  dev.off()
  
  #return the information
  coef = data.frame(yint,b1c,b2c,b3c,b4c)
  tbls = list(coef,final,g)
  return(tbls)
}
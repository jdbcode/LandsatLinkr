#' Predict oli TC index
#'
#' Predict oli TC index
#' @param file The full path name of the *archv file
#' @param pngout Logical. If TRUE then a 4-panel .png image file will be created that for each oli band displays a truncated histogram and the estimated Dark Object Value
#' @import MASS
#' @export


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
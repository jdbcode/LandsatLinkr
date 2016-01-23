#' Predict MSS TC index
#'
#' Predict MSS TC index
#' @param file The full path name of the *archv file
#' @param pngout Logical. If TRUE then a 4-panel .png image file will be created that for each MSS band displays a truncated histogram and the estimated Dark Object Value
#' @import MASS
#' @export


predict_mss_index = function(tbl, outsampfile){  
  #predict_mss_index = function(refsamp, b1samp, b2samp, b3samp, b4samp, mss_sr_file, ref_file, index, outsampfile, samplen){
  
  #tbl = data.frame(mss_img = rep(basename(mss_sr_file),samplen),
  #                 ref_img = rep(basename(ref_file),samplen),
  #                 index = rep(index,samplen),
  #                 refsamp,b1samp,b2samp,b3samp,b4samp)
  #final = tbl[complete.cases(tbl),]

  #create a multivariable linear model
  model = rlm(refsamp ~ b1samp + b2samp + b3samp + b4samp, data=tbl) #tbl replaced final 1/22/2016
  
  #final$singlepred = round(predict(model))
  #write.csv(final, outsampfile, row.names=F)
  
  tbl$singlepred = round(predict(model))
  write.csv(tbl, outsampfile, row.names=F)
  
  #apply the model to the original img
  #yint = model$coefficients[1]
  #b1c = model$coefficients[2]
  #b2c = model$coefficients[3]
  #b3c = model$coefficients[4]
  #b4c = model$coefficients[5]
  
  #plot the regression
  #r = cor(final$refsamp, final$singlepred)
  #coef = rlm(final$refsamp ~ final$singlepred)
  r = cor(tbl$refsamp, tbl$singlepred)
  coef = rlm(tbl$refsamp ~ tbl$singlepred)
  
  pngout = sub("samp.csv", "plot.png",outsampfile)
  png(pngout,width=700, height=700)
  title = paste(tbl$index[1],"linear regression: slope =",paste(signif(coef$coefficients[2], digits=3),",",sep=""),
                "y Intercept =",paste(round(coef$coefficients[1], digits=3),",",sep=""),
                "r =",signif(r, digits=3))
  plot(x=tbl$singlepred,y=tbl$refsamp, #tbl replaced final 1/22/2016
       main=title,
       xlab=paste(tbl$mss_img[1],tbl$index[1]), #basename(mss_sr_file)    tbl$index[1] replaced index 1/22/2016
       ylab=paste(tbl$ref_img[1],tbl$index[1]))    #basename(ref_file)     tbl$index[1] replaced index 1/22/2016
  abline(coef = coef$coefficients, col="red")  
  dev.off()
  
  #return the information
  #coef = data.frame(yint,b1c,b2c,b3c,b4c)
  #tbls = list(coef,final) #tbl replaced final 1/22/2016
  #return(tbls)
  
  coef_tbl = data.frame(rbind(model$coefficients))
  cnames = c("yint","b1c","b2c","b3c","b4c")
  colnames(coef_tbl) = cnames
  tbls = list(coef_tbl,tbl)
  return(tbls)
  
}
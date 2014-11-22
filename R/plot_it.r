#' scatterplot of MSS and TM TC transformations
#'
#' create scatterplot of calibrated MSS and TM TC transformations
#' @param samp data.frame dataframe containing sample of pixel values for MSS and TM TC transformations 
#' @param coef array linear regression coefficents for linear model between MSS and TM TC transformations
#' @param r correlation coefficent for linear relationship between MSS and TM TC transformations
#' @param mss_file filename full path to MSS file
#' @param ref_file filename full path to TC reference file
#' @param index character what index to plot
#' @import ggplot2

plot_it = function(samp, coef, r, mss_file, ref_file, index){
  g = ggplot(samp, aes(msssamp, refsamp)) +
    stat_binhex(bins = 100)+
    scale_fill_gradientn(name = "Count", colours = rainbow(7))+
    xlab(paste(basename(mss_file),index)) +
    ylab(paste(basename(ref_file),index)) +
    ggtitle(paste("Linear Regression: Slope =",paste(signif(coef$coef, digits=3),",",sep=""),
                  "Y Intercept =",paste(round(coef$yint, digits=3),",",sep=""),
                  "r =",signif(r, digits=3))) +
    theme(plot.title = element_text(size = 12)) +
    geom_smooth(method="rlm", colour = "black")
  return(g)
}
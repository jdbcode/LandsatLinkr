#' Create TC plane plots
#'
#' Create TC plane plots
#' @param df data.frame. data.frame with column names "x" and "y" to be displayed as a 2d scatterplot 
#' @param plane character. the TC plane to be plotted options: "bg", "bw", "wg"
#' @param sensor character. which sensor is this data from options: "mss", "tm"
#' @import ggplot2
#' @export


plot_tc_planes = function(df,plane, sensor){
  blimits = c(-500, 10000)
  glimits = c(-500, 5000)
  wlimits = c(-6000, 1000)
  if(plane == "bg" | plane == "bw"){
    xlim = blimits
    xlab = paste(sensor, "tc brightness")
  }
  if(plane == "wg"){
    xlim = wlimits
    xlab = paste(sensor, "tc wetness")
  }
  if(plane == "bg" | plane == "wg"){
    ylim = glimits
    ylab = paste(sensor, "tc greenness")
  }
  if(plane == "bw"){
    ylim = wlimits
    ylab = paste(sensor, "tc wetness") 
  }
  
  tcplot = ggplot(df, aes(x, y)) +
    stat_binhex(bins = 100)+
    scale_fill_gradientn(colours = rainbow(7)) +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = ylim) +
    xlab(xlab) +
    ylab(ylab) +
    theme(legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1))+
    theme_bw()
  
  return(tcplot)
}
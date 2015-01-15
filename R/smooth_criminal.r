
library(raster)
library(segmented)
library(zoo)


ls_profile_stack = function(){
  x = -1174943
  y =  2528658
  index = "tca"
  dir = "L:/composites/test9"

  indexdir = file.path(dir,index)
  img = list.files(indexdir, "stack.bsq$", full.name=T)
  year = as.numeric(unique(substr(basename(list.files(indexdir, "composite.bsq$")),1,4)))
  r = brick(img)
  
  if(index == "tcb"){limits = c(-500, 10000)}
  if(index == "tcg"){limits = c(-500, 5000)}
  if(index == "tcw"){limits = c(-6000, 1000)}
  if(index == "tca"){limits = c(-500, 5000)}
  
  x = -1181263
  y = 2537648
  xy = data.frame(x=x,y=y)
  origv = value = as.vector(extract(r,xy))
  
  before = c(NA, value[-1] - value[-length(value)])*-1
  after = c(value[-1] - value[-length(value)], NA)
  absba = rowSums(cbind(abs(before),abs(after)), na.rm=T)
  difba = abs(rowSums(cbind(before,after*-1), na.rm=T))
  
  ratio = difba/absba
  these = which(ratio <= 0.15)
  
  #medabsba = median(absba,na.rm=T)
  #maddifba = median(difba,na.rm=T)
  
  #these = which(absba > medabsba & difba < maddifba)
  adj = rowMeans(cbind(before[these],after[these]))
  value[these] = adj+value[these]
  
#   value[value==0] = NA
#   nas = is.na(value)
#   goods = which(nas == F)
#   if(nas[1] == T)e[1] = e[goods[1]]
#   if(nas[n] == T)e[n] = e[goods[length(goods)]]
#   e = na.approx(e)
  
  df = data.frame(value,year)
#   mylm = lm(value ~ year)
  mylm = lm(value ~ year, data = df)
  myseg = 7
  while (class(myseg[1]) == "numeric"){
    myseg = ls_segit(mylm, seg.Z = ~ year, psi = NA, seg.control(stop.if.error = F, K = myseg, n.boot=0,it.max=100)) #it.max=20
  }
  df$fitted = round(fitted(myseg))
  
g = ggplot() +
  geom_point(data=df, aes(x=year, y=value), size = 4) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_line(data=df, aes(x=year, y=fitted)) +
  ylim(limits) +
  scale_x_continuous(breaks = seq(1972,2014))+ 
  ylab(index)

y1 = e.divisive(X=matrix(value),sig.lvl=0.1,R=199,k=NULL,min.size=2,alpha=1)
estimate = y1$estimate
estimate[length(estimate)] = estimate[length(estimate)]-1
fitted = data.frame(year = year[estimate], value = value[estimate])

g = ggplot() +
  geom_point(data=df, aes(x=year, y=value), colour = "red") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_point(data=fitted, aes(x=year, y=value, size = 6)) +
  ylim(limits) +
  scale_x_continuous(breaks = seq(1972,2014))+ 
  ylab(index)



y2 = e.agglo(X=matrix(value))

library(ecp)  
  
  

  
  
  
  
  
}
  

index = "tca"
dir = "L:/composites/test9"
indexdir = file.path(dir,index)
img = list.files(indexdir, "stack.bsq$", full.name=T)
r = brick(img)
year = as.numeric(unique(substr(basename(list.files(indexdir, "composite.bsq$")),1,4)))

if(index == "tcb"){limits = c(-500, 10000)}
if(index == "tcg"){limits = c(-500, 5000)}
if(index == "tcw"){limits = c(-6000, 1000)}
if(index == "tca"){limits = c(-500, 5000)}


x = -1163688
y =  2522203
year = as.numeric(unique(substr(basename(list.files(indexdir, "composite.bsq$")),1,4)))

smooth_criminal = function(x,y,r,year){
  value = as.vector(sampleRandom(r, size=1, na.rm=TRUE, xy=T))
  value = value[3:length(value)]
  #xy = data.frame(x=x,y=y)
#value = as.vector(extract(r,xy))  

#despike
for(i in 1:2){
before = c(NA, value[-1] - value[-length(value)])*-1
after = c(value[-1] - value[-length(value)], NA)
absba = rowSums(cbind(abs(before),abs(after)), na.rm=T)
difba = abs(rowSums(cbind(before,after*-1), na.rm=T))
ratio1 = absba/median(absba)
ratio2 = difba/median(difba)
these = which(ratio1 > 4 & ratio2 < 5)
adj = rowMeans(cbind(before[these],after[these]))
value[these] = adj+value[these]
}
df = data.frame(value,year,fitted=NA)

#breakpoints
y1 = e.divisive(X=matrix(value),sig.lvl=0.05,R=199,k=NULL,min.size=2,alpha=1)
estimate = y1$estimate
estimate[length(estimate)] = estimate[length(estimate)]-1

#add fast bp
for(i in 1:3){
  befm = median(before[before > 0], na.rm=T)
  #aftm = median(after[after < 0], na.rm=T)
  these = which(before[estimate] > befm)
  estimate = c(estimate, estimate[these]-1)
}
estimate = sort(unique(estimate))
# fitted = data.frame(year = year[estimate], value = value[estimate])

# #plot
# g = ggplot() +
#   geom_point(data=df, aes(x=year, y=value), colour = "red") +
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
#   geom_point(data=fitted, aes(x=year, y=value, size = 6)) +
#   ylim(limits) +
#   scale_x_continuous(breaks = seq(1972,2014))+ 
#   ylab(index)


len = length(estimate)
breaks = array(NA, len)
for(i in 1:(len-1)){
  print(i)
  these = seq(estimate[i], estimate[i+1])
  newdf = df[these,]
  line = lm(value ~ year, newdf)
  v = fitted(line)
  breaks[i] = mean(c(breaks[i], v[1]),na.rm=T)
  breaks[i+1] = v[length(these)]
  #df$fitted[these] = round(fitted(line))
}
for(i in 1:(len-1)){
  print(i)
  value = breaks[c(i,i+1)]
  year = df$year[c(estimate[i],estimate[i+1])]
  these = seq(estimate[i], estimate[i+1])
  line = lm(value ~ year)
  v = predict(line, df[these,])
  df$fitted[these] = v
}


g = ggplot() +
  geom_point(data=df, aes(x=year, y=value), colour = "red") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_line(data=df, aes(x=year, y=fitted)) +
  ylim(limits) +
  scale_x_continuous(breaks = seq(1972,2014))+ 
  ylab(index)

print(g)

}
  
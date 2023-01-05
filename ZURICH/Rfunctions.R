my.germany.map = function(data, cutpoints=seq(min(data),max(data),length=256), autoscale=FALSE, legend=TRUE, append=FALSE)
{
  if (autoscale)
  {
    data = (data-min(data))/(max(data)-min(data)+1e-8)
  }
  #cutpoints = c(-1e9,cutpoints, 1e9)
  
  # farben <- rainbow(as.numeric(cut(data,cutpoints,include.lowest=T))/length(cutpoints))
  cols =  hcl.colors(256)
  farben <- cols[as.numeric(cut(data,cutpoints,include.lowest=T))]
  
  xmin <- 1:length(germany)
  xmax <- 1:length(germany)
  ymin <- 1:length(germany)
  ymax <- 1:length(germany)
  
  for(i in 1:length(germany))
  {
    xmin[i] <- min(germany[[i]][,2],na.rm=T)
    xmax[i] <- max(germany[[i]][,2],na.rm=T)
    ymin[i] <- min(germany[[i]][,3],na.rm=T)
    ymax[i] <- max(germany[[i]][,3],na.rm=T)
  }
  
  breite <- c(min(xmin),max(xmax))
  hoehe <- c(min(ymin),max(ymax))
  
  if (!append) plot(breite,hoehe,type="n",axes=F, xlab=" ", ylab=" ", asp = 1)
  
  
  for(k in length(germany):1)
  {
    polygon(germany[[k]][,2],germany[[k]][,3],col=farben[k])
  }
  
  
}
library(tidyverse)
library(INLA)
library(inlabru)
library(ggplot2)
data(meuse)
data(meuse.grid)
meuse  = meuse %>% dplyr::select(x,y,zinc, dist)


mesh = inla.mesh.2d(loc.domain = cbind(meuse$x, meuse$y), 
                     max.edge = c(150, 500),
                     cutoff = 100,
                     offset = c(100, 1250) )


meuse1 = meuse %>% dplyr::select(x,y,dist)
coordinates(meuse1) = c("x","y")
grd <- expand.grid(x=seq(from=min(meuse$x)-200, to=max(meuse$x)+200, length.out = 200), 
                   y=seq(from=min(meuse$y)-200, to=max(meuse$y)+200, length.out = 200)) 
coordinates(grd) <- ~ x+y
gridded(grd) <- TRUE

template = raster::rasterFromXYZ(grd)

library(fields)
spline = Tps(meuse1@coords, meuse1$dist)
splined = raster::interpolate(template, spline)
plot(splined)

aa  = data.frame(x = coordinates(splined)[,1],
                 y = coordinates(splined)[,2],
                 z = raster::values(splined))

aa$z[aa$z<0] = -aa$z[aa$z<0]

ggplot() + geom_tile(data = aa, aes(x,y,fill = z)) +
  gg(mesh)

xx = range(mesh$loc[,1])
yy = range(mesh$loc[,2])

new_xy  = as.matrix(expand.grid(seq(xx[1],xx[2],length.out =  400),
                                seq(yy[1],yy[2],length.out =  400)))

A = inla.spde.make.A(mesh = mesh,loc = new_xy)
ips <- ipoints(domain = mesh)


sst_SPDF = SpatialPixelsDataFrame(aa[,c(1,2)], 
                                  data = data.frame(z =aa[,3]))

v <- inlabru:::eval_SpatialDF(data = sst_SPDF,
                              where = ips)
v2 <- inlabru:::bru_fill_missing(data = sst_SPDF,
                                 where = ips,
                                 values = v)


val   = as.vector(A%*% v2)

data.frame(x = new_xy[,1], 
           y = new_xy[,2], 
           dist = val) %>%
  ggplot() + geom_tile(aes(x,y,fill = val))


data(meuse.grid)
coordinates(meuse.grid) = ~x+y
gridded(meuse.grid) = TRUE

boundary = maptools::unionSpatialPolygons(
  as(meuse.grid, "SpatialPolygons"), rep (1, length(meuse.grid))
)

meuseData = list(meuse  = meuse,
                 dist_raster = data.frame(x = new_xy[,1], 
                                          y = new_xy[,2], 
                                          dist = val),
                 boundary = boundary)
saveRDS(meuseData, "Data/meuseData")

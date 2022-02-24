library(tidyverse)
library(ggplot2)
library(inlabru)
library(INLA)
Piemonte_data = read.table("Day2.1/Data/Piemonte_data/Piemonte_data_byday.csv",
                           sep = ",", header = T)

nstat = length(unique(Piemonte_data$UTMX))
nday = length(unique(Piemonte_data$Date))


Piemonte_data$Date = as.Date(Piemonte_data$Date, format = "%d/%m/%y" )
Piemonte_data$time = rep(1:nday, each = nstat)



df = Piemonte_data %>% 
  rename(x = UTMX,
         y = UTMY,
         dem =  A) %>%
  mutate(logPM10 = log(PM10)) %>%
  dplyr::filter(time<20)


border = read.table("Day2.1/Data/Piemonte_data/Piemonte_borders.csv",
                   header = T, sep = ",")
mesh = inla.mesh.2d(loc = cbind(df$x, df$y),
                    loc.domain = border,
                    offset = c(20,40),
                    max.edge = c(30,50))
plot(mesh)
lines(border)
spde = inla.spde2.pcmatern(mesh = mesh,
                           prior.range=c(100,0.5), 
                           prior.sigma=c(1,0.1)) 



coordinates(df) = c("x","y")
cmp  = ~ Intercept(1) + 
  SPDE(coordinates, model = spde,
       group = time, control.group = list(model = "ar1")) +
  alt(dem, model = "linear") +
  random(Station.ID , model = "iid")

lik = like(formula = logPM10 ~ Intercept + SPDE + alt,
           family = "gaussian",
           data = df)

fit = bru(cmp, lik,
          options = list(verbose = F,
                         bru_max_iter = 1,
                         inla.mode  = "experimental"))

library(sf)
bb = sf::st_polygon(list(as.matrix(border))) %>%
  st_simplify(dTolerance = 0.1)

shape_sp = as(bb,"Spatial")
pxl = pixels(mesh, nx = 200, ny = 200, mask = shape_sp)
ips2 <- ipoints(domain = c(1:19), name = "time")

pxl_time = cprod(pxl,ips2)
pred = predict(fit, pxl_time)


load("Day2.1/Data/Piemonte_data/Covariates/Altitude_GRID.Rdata")
load("Day2.1/Data/Piemonte_data/Covariates/Piemonte_grid.Rdata")


data.frame(x = Piemonte_grid[,1],
           y = Piemonte_grid[,2],
           z = as.vector(t(AltitudeGRID))) %>%
  ggplot() + geom_tile(aes(x,y,fill = z))


## object with altiture, this does not cover the whole mesh!
dem = SpatialPixelsDataFrame(points = cbind(Piemonte_grid[,1],
                                            Piemonte_grid[,2]),
                             data = data.frame(dem = as.vector(t(AltitudeGRID))))

saveRDS(dem, file = "dem.RDS")



## original data, this does not cover the mesh
dem = readRDS("dem.RDS")
## create another object which  does cover the whole mesh
large_grid = expand_grid(x = seq(250, 580,4),
                         y = seq(4810, 5210,4))
dem_large = SpatialPixelsDataFrame(points = large_grid,
                                   data = data.frame(dem = 
                                                       rep(NA,length(large_grid$x))))

# Use bru_fill_missing to fill the new object
v2 <- bru_fill_missing(data = dem,
                       where = dem_large,
                       values = dem$dem)


dem_large$dem = v2


bru_fill_missing(dem, pxl)

predict(fit, pxl, 
        alt_eval(inlabru:::eval_SpatialDF(covar_grid, .data.))



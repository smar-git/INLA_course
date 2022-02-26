library(tidyverse)
library(ggplot2)
library(inlabru)
library(INLA)
library(viridis)

## Load the data
df = read.table("Day2.1/Data/Piemonte_data/Piemonte_data.dat")

df = df[df$time<=10,]

# area of interest
shape = readRDS("Day2.1/Data/Piemonte_data/Piemonte")

# create the mesh and hte spde model
mesh = inla.mesh.2d(loc = cbind(df$x, df$y),
                    offset = c(20,40),
                    max.edge = c(30,50))

spde = inla.spde2.pcmatern(mesh = mesh,
                           prior.range=c(100,0.5), 
                           prior.sigma=c(1,0.1)) 


# make the data a spatial object
coordinates(df) = c("x","y")

# model component
cmp  = ~ Intercept(1) + 
  SPDE(coordinates, model = spde,
       group = time, control.group = list(model = "ar1")) +
  dem(dem, model = "linear")

# likelihood
lik = like(formula = logPM10 ~ Intercept + SPDE + dem + temp ,
           family = "gaussian",
           data = df)

# fit the model
fit = bru(cmp, lik,
          options = list(verbose = F,
                         bru_max_iter = 1,
                         inla.mode  = "experimental"))

fit$summary.fixed
## predictions

pxl = pixels(mesh, nx = 200, ny = 200, mask = shape)
ips2 <- ipoints(domain = c(1:5), name = "time")

pxl_time = cprod(pxl,ips2)

alt = read.table("Day2.1/Data/Piemonte_data/Altitude.dat")

## object with altiture, this does not cover the whole mesh!
dem = SpatialPixelsDataFrame(points = cbind(alt[,1],
                                            alt[,2]),
                             data = data.frame(dem = alt[,3]))


## create another object which  does cover the whole mesh
large_grid = expand_grid(x = seq(250, 580,4),
                         y = seq(4810, 5210,4))


# expand dem --------------------------------------------------------------
dem_large = SpatialPixelsDataFrame(points = large_grid,
                                   data = data.frame(dem = 
                                                       rep(NA,length(large_grid$x))))

# Use bru_fill_missing to fill the new object
dem_large$dem <- bru_fill_missing(data = dem,
                       where = dem_large,
                       values = dem_large$dem)



# expand temp -------------------------------------------------------------
 load("Day2.1/Data/Piemonte_data/Temp_GRID.Rdata")
temp = Mean_Temp[,,1:5]

mat = matrix(NA, dim(temp)[1]*dim(temp)[2], 5)
for(i in 1:5)
  mat[,i] =  as.vector(t(temp[,,i]))
data.frame(x =alt[,1],
           y =alt[,2],
           
           temp = mat) %>%
  pivot_longer(-c(x,y)) %>%
  ggplot() + geom_tile(aes(x,y,fill = value)) +
  facet_wrap(.~name)


temp = SpatialPixelsDataFrame(points = alt[,c(1,2)],
                             data = data.frame(temp = mat))

temp_large = SpatialPixelsDataFrame(points = cbind(large_grid$x,
                                                   large_grid$y),
                                   data = data.frame(temp = matrix(NA, 
                                                                   length(large_grid$x),5)))
                                                       


# # Use bru_fill_missing to fill the new objectt

fill_names = names(temp)
for (layer in fill_names) {
  print(layer)
  temp_large[[layer]] <- bru_fill_missing(temp, 
                                          temp_large, 
                                          values = temp_large[[layer]], 
                                          layer = layer)
}



# prediction at stations--------------------------------------------------------------

 
pred_at_station = predict(fit, df, ~exp( Intercept + 
                        SPDE + 
                        dem ))

sel = sample(1:24, 6) %>% sort()
as.data.frame(pred_at_station) %>%
  dplyr::filter(Station.ID %in% sel) %>% ggplot() + 
  geom_line(aes(time, PM10, group = Station.ID ), color = "red") +
  geom_line(aes(time, median, group = Station.ID)) +
  geom_ribbon(aes(time, 
                  ymin= q0.025, 
                  ymax = q0.975,
                  group = Station.ID), alpha = 0.5) +
  facet_wrap(.~Station.ID)


# prediction in space -----------------------------------------------------

pxl_time <- cbind(cprod(pxl,ips2), temp = 0)

pred_space = generate(fit, pxl_time, 
                     ~temp_eval(inlabru:::eval_SpatialDF(temp_large,
                                                           .data.,
                                                         selector = "time")),
                     n.samples = 1)


pred_space = predict(fit, pxl_time, 
        ~data.frame( 
          spde = SPDE,
          log_scale = Intercept + 
                       SPDE + 
                       dem_eval(inlabru:::eval_SpatialDF(dem_large,
                                                         .data.)),
          nat_scale = exp( Intercept + 
                             SPDE + 
                             dem_eval(inlabru:::eval_SpatialDF(dem_large, 
                                                               .data.)))))

data.frame(x = coordinates(pxl_time)[,1],
           y = coordinates(pxl_time)[,2],
           time = pxl_time$time,
           mean = pred_space$log_scale$sd) %>%
  ggplot() + geom_tile(aes(x,y,fill  = mean)) +
  scale_fill_viridis() +
  facet_wrap(.~time)

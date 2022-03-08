## ----sett, echo = FALSE, results = 'hide', message = FALSE, warning = FALSE----
source("R/initial_setup.R")
opts_chunk$set(
  fig.path = 'figs/prefsampl-'
)
library(splancs)

## ----window, message = FALSE, warning = FALSE----------------------------
library(spatstat)
win <- owin(c(0, 3), c(0, 3))

## ----gridres-------------------------------------------------------------
npix <- 300
spatstat.options(npixel = npix)

## ----n-exp---------------------------------------------------------------
beta0 <- 3
exp(beta0) * diff(range(win$x)) * diff(range(win$y))

## ------------------------------------------------------------------------
sigma2x <- 0.2
range <- 1.2
nu <- 1

## ----simulapp,eval=TRUE, warning=FALSE, message=FALSE--------------------
library(RandomFields)
set.seed(1)
lg.s <- rLGCP('matern', beta0, var = sigma2x,
  scale = range / sqrt(8), nu = nu, win = win)

## ----xy------------------------------------------------------------------
xy <- cbind(lg.s$x, lg.s$y)[, 2:1]

## ----nxy-----------------------------------------------------------------
(n <- nrow(xy))

## ------------------------------------------------------------------------
Lam <- attr(lg.s, 'Lambda')
rf.s <- log(Lam$v)
summary(as.vector(rf.s))

## ----lgrfpp, echo = FALSE, fig.cap = "Simulated intensity of the point process and simulated point pattern (black dots)."----
par(mfrow =c(1, 1), mar=c(0, 0, 0, 0))
book.plot.field(list(x = Lam$yrow, y = Lam$xcol, z = rf.s), 
                xlim=c(0,3), ylim=c(0,3))
points(xy, pch = 19) 

## ----mesh----------------------------------------------------------------
loc.d <- 3 * cbind(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0))
mesh <- inla.mesh.2d(loc.domain = loc.d, offset = c(0.3, 1), 
  max.edge = c(0.3, 0.7), cutoff = 0.05)
nv <- mesh$n

## ----ppmesh, echo = FALSE, fig.cap = "Mesh used to fit a log-Gaussian Cox process to a point pattern."----

par(mar = c(0, 0, 0, 0))
plot(mesh, asp = 1, main = '')
points(xy, col = 4, pch = 19)
lines(loc.d, col = 3)

## ----spde----------------------------------------------------------------
spde <- inla.spde2.pcmatern(mesh = mesh,
  # PC-prior on range: P(practic.range < 0.05) = 0.01
  prior.range = c(0.05, 0.01),
  # PC-prior on sigma: P(sigma > 1) = 0.01
  prior.sigma = c(1, 0.01)) 

## ----dualmesh0-----------------------------------------------------------
dmesh <- book.mesh.dual(mesh)

## ----splocd--------------------------------------------------------------
domain.polys <- Polygons(list(Polygon(loc.d)), '0')
domainSP <- SpatialPolygons(list(domain.polys))

## ----pols, message=FALSE, warning=FALSE----------------------------------
library(rgeos)
w <- sapply(1:length(dmesh), function(i) {
  if (gIntersects(dmesh[i, ], domainSP))
    return(gArea(gIntersection(dmesh[i, ], domainSP)))
  else return(0)
})

## ------------------------------------------------------------------------
sum(w)

## ----wsummary------------------------------------------------------------
table(w > 0)

## ----dualmesh, echo = FALSE, fig.cap = "Voronoy polygons for the mesh used to make inference for the log-Gaussian Cox process."----

par(mar = c(2, 2, 1, 1), mgp = 2:0)
plot(mesh$loc, asp = 1, col = (w == 0) + 1, pch = 19, xlab = '', ylab = '') 
plot(dmesh, add = TRUE)
lines(loc.d, col = 3)

## ----y01-----------------------------------------------------------------
y.pp <- rep(0:1, c(nv, n))

## ----expected------------------------------------------------------------
e.pp <- c(w, rep(0, n)) 

## ----pp-proj-------------------------------------------------------------
imat <- Diagonal(nv, rep(1, nv))

## ----Aloc----------------------------------------------------------------
lmat <- inla.spde.make.A(mesh, xy)

## ----App-----------------------------------------------------------------
A.pp <- rbind(imat, lmat)

## ----stkpp---------------------------------------------------------------
stk.pp <- inla.stack(
  data = list(y = y.pp, e = e.pp), 
  A = list(1, A.pp),
  effects = list(list(b0 = rep(1, nv + n)), list(i = 1:nv)),
  tag = 'pp')

## ----ppest---------------------------------------------------------------
pp.res <- inla(y ~ 0 + b0 + f(i, model = spde), 
  family = 'poisson', data = inla.stack.data(stk.pp), 
  control.predictor = list(A = inla.stack.A(stk.pp)), 
  E = inla.stack.data(stk.pp)$e)

## ----pppars,eval=TRUE, R.options = list(digits = 3)----------------------
pp.res$summary.hyperpar

## ----pppost, echo = FALSE, results = 'hide', fig.width = 8, fig.height = 4, fig.cap = '(ref:pppost)'----

par(mfrow = c(1, 3), mar = c(3, 3, 1, 0.3), mgp = c(2, 1, 0)) 

plot(pp.res$marginals.fixed[[1]], type = 'l', xlab = expression(beta[0]),
     ylab = 'Density')
abline(v = beta0, col = 2)

plot(pp.res$marginals.hyperpar[[2]], type = 'l', xlab = expression(sigma),
     ylab = 'Density', xlim = c(0,2))
abline(v = sqrt(sigma2x), col = 2)

plot(pp.res$marginals.hyperpar[[1]], type = 'l', xlab = 'Nominal range',
     ylab = 'Density', xlim = c(0, 8))
abline(v = range, col = 2)

## ----gridcov-------------------------------------------------------------
# Use expanded range
x0 <- seq(min(mesh$loc[, 1]), max(mesh$loc[, 1]), length = npix)
y0 <- seq(min(mesh$loc[, 2]), max(mesh$loc[, 2]), length = npix)
gridcov <- outer(x0, y0, function(x,y) cos(x) - sin(y - 2))

## ----n-exp-cov-----------------------------------------------------------
beta1 <- -0.5
sum(exp(beta0 + beta1 * gridcov) * diff(x0[1:2]) * diff(y0[1:2]))

## ----simulappc-----------------------------------------------------------
set.seed(1)
lg.s.c <- rLGCP('matern', im(beta0 + beta1 * gridcov, xcol = x0,
  yrow = y0), var = sigma2x, scale = range / sqrt(8), 
  nu = 1, win = win)

## ----xyc-----------------------------------------------------------------
xy.c <- cbind(lg.s.c$x, lg.s.c$y)[, 2:1]
n.c <- nrow(xy.c)

## ----lgrfppc, echo = FALSE, results = 'hide', fig.width = 9, fig.height = 4, fig.cap = "Simulated covariate (left) and simulated log-intensity of the point process, along with the simulated point pattern (right)."----
par(mfrow = c(1, 2), mar = c(2, 1, 0.5, 4.5), mgp = c(1, 0.5, 0))
book.plot.field(list(x = x0, y = y0, z = gridcov), xlim = c(0, 3),
  ylim = c(0, 3))
book.plot.field(list(x = x0, y = y0, z = log(attr(lg.s.c, 'Lambda')$v)),
                xlim = c(0, 3), ylim = c(0, 3))  
points(xy.c, pch = 19)

## ----collcovar, results = "hide"-----------------------------------------
covariate.im <- im(gridcov, x0, y0)
covariate <- interp.im(covariate.im, 
  x = c(mesh$loc[, 1], xy.c[, 1]),
  y = c(mesh$loc[, 2], xy.c[, 2]))

## ----datc----------------------------------------------------------------
y.pp.c <- rep(0:1, c(nv, n.c))
e.pp.c <- c(w, rep(0, n.c))

## ----A.c-----------------------------------------------------------------
lmat.c <- inla.spde.make.A(mesh, xy.c)

## ----App.c---------------------------------------------------------------
A.pp.c <- rbind(imat, lmat.c)

## ----stkpp.c-------------------------------------------------------------
stk.pp.c <- inla.stack(
  data = list(y = y.pp.c, e = e.pp.c), 
  A = list(1, A.pp.c), 
  effects = list(list(b0 = 1, covariate = covariate), 
    list(i = 1:nv)),
  tag = 'pp.c')

## ----ppest.c-------------------------------------------------------------
pp.c.res <- inla(y ~ 0 + b0 + covariate + f(i, model = spde),
  family = 'poisson', data = inla.stack.data(stk.pp.c), 
  control.predictor = list(A = inla.stack.A(stk.pp.c)), 
  E = inla.stack.data(stk.pp.c)$e)

## ----insp-u.c, R.options = list(digits = 3)------------------------------
pp.c.res$summary.hyperpar

## ----pppostc, echo = FALSE, results = 'hide', fig.cap = '(ref:pppostc)'----

par(mfrow = c(2, 2), mar = c(3, 3, 1, 0.3), mgp = c(2, 1, 0)) 

plot(pp.c.res$marginals.fix[[1]], type = 'l', ylab = 'Density', 
  xlab = expression(beta[0]))
abline(v = beta0, col = 2)

plot(pp.c.res$marginals.fix[[2]], type = 'l', ylab = 'Density', 
  xlab = expression(beta[1]))
abline(v = beta1, col = 2)

plot(pp.c.res$marginals.hyperpar[[2]], type = 'l', ylab = 'Density', 
  xlab = expression(sigma), xlim = c(0, 2))
abline(v = sqrt(sigma2x), col = 2)

plot(pp.c.res$marginals.hyperpar[[1]], type = 'l', ylab = 'Density',
  xlab = "Spatial range", xlim = c(0, 10))
abline(v = range, col = 2)





## ----simulaz-------------------------------------------------------------
z <- log(t(Lam$v)[do.call('cbind',
  nearest.pixel(xy[, 1], xy[, 2], Lam))])

## ------------------------------------------------------------------------
summary(z)

## ----resp----------------------------------------------------------------
beta0.y <- 10
beta <- -2
prec.y <- 16

set.seed(2)
resp <- beta0.y + (z - beta0) / beta + 
  rnorm(length(z), 0, sqrt(1 / prec.y))

## ------------------------------------------------------------------------
summary(resp)

## ----rresp---------------------------------------------------------------
stk.u <- inla.stack(
  data = list(y = resp),
  A = list(lmat, 1), 
  effects = list(i = 1:nv, b0 = rep(1, length(resp))))

u.res <- inla(y ~ 0 + b0 + f(i, model = spde),
  data = inla.stack.data(stk.u), 
  control.predictor = list(A = inla.stack.A(stk.u)))


# inlabru() ---------------------------------------------------------------

df_points = data.frame(x = xy[,1],
                       y = xy[,2])
coordinates(df_points) = c("x","y")
df_marks = data.frame(x = xy[,1],
           y = xy[,2],
           z = resp)
coordinates(df_marks) = c("x","y")

cmp = ~Intercept1(1) + Intercept2(1) +
  field(coordinates, model = spde) +
  field_copy(coordinates, copy = "field", fixed = FALSE)

lik1 = like(formula = coordinates ~ Intercept1 + field,
            family = "cp",
            samplers = domainSP,
            domain = list(coordinates = mesh),
            data = df_points)

lik2 = like(formula = z ~ Intercept2  + field_copy,
            family = "gaussian",
            data = df_marks)

res = bru(cmp,
          lik1, lik2,
          options = list(verbose = F,
                         bru_verbose = 1,
                         inla.mode  = "experimental"))

## ----label = "ppres", echo = FALSE---------------------------------------
tab.ppres <- cbind(
  Parameters = c("$\\beta_0^y$", "$1/\\sigma^2_y$"), 
  True = c(beta0y = beta0.y, prec.y = prec.y), 
  rbind(u.res$summary.fixed[, c(1:3,5)], 
        u.res$summary.hyperpar[1, c(1:3,5)]))

#Column names
names(tab.ppres) <- c("Parameter", "True", "Mean", "St. Dev.",
  "2.5\\% quant.", "97.5\\% quant.")

knitr::kable(tab.ppres,
row.names = FALSE,
  caption = "Posterior modes of some of the model parameters.",
  format = "pandoc")

## ----upost, echo = FALSE, results = 'hide', fig.width = 6, fig.height = 2.1, fig.cap = '(ref:upost)'----

par(mfrow = c(1, 3), mar = c(3, 3, 0.3, 0.3), mgp = c(2, 1, 0))

marg.sigma <- inla.tmarginal(function(x) sqrt(1 / x),
  u.res$marginals.hyperpar[[1]])

plot(marg.sigma, type = 'l', ylab = 'Density', xlab = expression(sigma)) 
abline(v = sqrt(1 / prec.y), col = 2)

plot(u.res$marginals.hyperpar[[2]], type = 'l', ylab = 'Density', 
  xlab = 'Spatial range', xlim = c(0, 10))
abline(v = range, col = 2)

plot(u.res$marginals.hyperpar[[3]], type = 'l', 
  xlab = expression(sigma[x]), ylab = 'Density')
abline(v = sqrt(sigma2x), col = 2)

## ----ppstk---------------------------------------------------------------
stk2.y <- inla.stack(
  data = list(y = cbind(resp, NA), e = rep(0, n)), 
  A = list(lmat, 1),
  effects = list(i = 1:nv, b0.y = rep(1, n)),
  tag = 'resp2')

stk2.pp <- inla.stack(data = list(y = cbind(NA, y.pp), e = e.pp), 
  A = list(A.pp, 1),
  effects = list(j = 1:nv, b0.pp = rep(1, nv + n)),
  tag = 'pp2')

j.stk <- inla.stack(stk2.y, stk2.pp)

## ----j-res---------------------------------------------------------------
# Gaussian prior
gaus.prior <- list(prior = 'gaussian', param = c(0, 2))
# Model formula
jform <- y ~ 0 + b0.pp + b0.y + f(i, model = spde) +
  f(j, copy = 'i', fixed = FALSE,
    hyper = list(theta = gaus.prior))
# Fit model
j.res <- inla(jform, family = c('gaussian', 'poisson'), 
  data = inla.stack.data(j.stk),
  E = inla.stack.data(j.stk)$e,
  control.predictor = list(A = inla.stack.A(j.stk)))

## ----label = "ppres2", echo = FALSE--------------------------------------
tab.ppres2 <- cbind(
  Parameters = c("$\\beta_0$", "$\\beta_0^y$"),
  True = c(beta0 = beta0, beta0y = beta0.y),
  j.res$summary.fixed[, c(1:3, 5)])

#Column names
names(tab.ppres2) <- c("Parameter", "True", "Mean", "St. Dev.",
  "2.5\\% quant.", "97.5\\% quant.")

knitr::kable(tab.ppres2,
row.names = FALSE,
  caption = "Posterior modes of some of the model parameters under preferential sampling.",
  format = "pandoc")

## ----jpars, echo = FALSE, results = 'hide', fig.width = 5.5, fig.height = 4, fig.cap = '(ref:jpars)'----

par(mfrow = c(2, 3), mar = c(3, 3, 0.5, 0.5), mgp = c(2, 0.5, 0))

plot(j.res$marginals.fix[[1]], type = 'l', ylab = 'Density', 
  xlab = expression(beta[0]), lwd = 2) 
lines(pp.res$marginals.fix[[1]], lty = 2, lwd = 2)
abline(v = beta0, col = 2)

plot(j.res$marginals.fix[[2]], type = 'l', ylab = 'Density', 
  xlab = expression(beta[0]^y), lwd = 2) 
lines(u.res$marginals.fix[[1]], lty = 3, lwd = 5)
abline(v = beta0.y, col = 2)

plot(inla.tmarginal(function(x) 1 / exp(x), 
    j.res$internal.marginals.hyperpar[[1]]), lwd = 2, 
  type = 'l', ylab = 'Density', xlab = expression(sigma[y]^2)) 
lines(inla.tmarginal(function(x) 1 / exp(x), 
    u.res$internal.marginals.hy[[1]]),
  lwd = 5, lty = 3)
abline(v = 1 / prec.y, col = 2)

plot(j.res$marginals.hyperpar[[2]], type = 'l', xlim = c(0, 10),
  xlab = 'Spatial range', ylab = 'Density', lwd = 2)
lines(pp.res$marginals.hyperpar[[1]], lty = 2, lwd = 2)
lines(u.res$marginals.hyperpar[[2]], lty = 3, lwd = 5)
abline(v = range, col = 2)

plot(j.res$marginals.hyperpar[[3]], type = 'l', lwd = 2, xlim = c(0, 1),
     xlab = expression(sqrt(sigma[x]^2)), ylab = 'Density')
lines(pp.res$marginals.hyperpar[[2]], lty = 2, lwd = 2)
lines(u.res$marginals.hyperpar[[3]], lty = 3, lwd = 5)
abline(v = sigma2x^0.5, col = 2)

plot(j.res$marginals.hyperpar[[4]], type = 'l', xlab = expression(beta),
  ylab = 'Density', lwd = 2)
abline(v = beta, col = 2)
legend('topright', c('True value', 'Joint', 'Only PP', 'Only Y'), 
  col = c(2, 1, 1, 1), lty = c(1, 1, 2, 3), lwd = 0.5 * c(2, 2, 3, 5),
  bty = 'n', cex = 0.65)


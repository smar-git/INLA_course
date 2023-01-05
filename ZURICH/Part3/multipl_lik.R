library(tidyverse)
library(INLA)
library(inlabru)



# Several likelihood

## Example: Coregionalization model

# $$
#   \begin{aligned}
# y_1(s) & =\alpha_1 + u(s) +e_1(s)\\
# y_2(s) & =\alpha_2+\lambda u(s) + e_2(s)
# \end{aligned}
# $$
#   
#   where the $\alpha_k$  are intercepts, $u(s)$ is the  spatial effect, $\lambda$ is a weights for  spatial effects and $e_k(s)$ are uncorrelated error terms, with $k=1,2,3$.



# Intercept on reparametrized model
alpha <- c(-5, 3) 
# Random field marginal variances
m.var <- c(0.5, 0.4) 
# GRF range parameters:
range <- c(4)
# Copy parameters: reparameterization of coregionalization 
# parameters
beta <- c(0.7) 
# Standard deviations of error terms
e.sd <- c(0.3, 0.2)


n1 <- 99
n2 <- 100


loc1 <- cbind(runif(n1) * 10, runif(n1) * 5) 
loc2 <- cbind(runif(n2) * 10, runif(n2) * 5) 

book.rMatern <- function(n, coords, sigma=1, range, kappa = sqrt(8*nu)/range, variance = sigma^2, nu=1) {
  m <- as.matrix(dist(coords))
  m <- exp((1-nu)*log(2) + nu*log(kappa*m)-
             lgamma(nu))*besselK(m*kappa, nu)
  diag(m) <- 1
  return(drop(crossprod(chol(variance*m),
                        matrix(rnorm(nrow(coords)*n), ncol=n))))
}
set.seed(05101980)
u <- book.rMatern(1, rbind(loc1, loc2), range = range[1],
                  sigma = sqrt(m.var[1]))


set.seed(08011952)

y1 <- alpha[1] + u[1:n1] + rnorm(n1, 0, e.sd[1])
y2 <- alpha[2] + beta * u[n1 + 1:n2] + 
  rnorm(n2, 0, e.sd[2])

mesh <- inla.mesh.2d(rbind(loc1, loc2), 
                     max.edge = c(0.5, 1.25), offset = c(0.1, 1.5), cutoff = 0.1)

plot(mesh)

spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(0.5, 0.01), # P(range < 0.5) = 0.01
  prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01

data1 = data.frame(x = loc1[,1],
                   y = loc1[,2],
                   z = y1)
coordinates(data1) = c("x","y")

data2 = data.frame(x = loc2[,1],
                   y = loc2[,2],
                   z = y2)
coordinates(data2) = c("x","y")

cmp = ~  Intercept1(1) + Intercept2(1) +
  u1(coordinates, model = spde) +
  u1_copy(coordinates, copy = "u1", fixed = FALSE)

lik1 = like(formula = z ~ Intercept1 + u1,
            family  = "gaussian",
            data = data1)

lik2 = like(formula = z ~ Intercept2 + u1_copy,
            family = "gaussian",
            data = data2)


res = bru(cmp, lik1, lik2)



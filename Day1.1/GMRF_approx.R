xi <- seq(-1.5,4,0.001)
y <- 3
delta <- 0.005
eta <- 0
approx <- function(xi, y, eta, delta){
  res <- dnorm(xi,eta, sqrt(1/delta))*dpois(y, exp(xi))
  return(res)
}
nominator_approx <- integrate(approx, -Inf, Inf, y=y, eta=eta,delta=delta)$value

# full conditional
fc <- function(xi, y, eta, delta){
  res <- exp(-(delta/2) * (xi - eta)^2 + y*xi - exp(xi))
  norm <-sqrt(delta)/(sqrt(2*pi) *factorial(y))
  #norm <- 1
  return((res *norm))
} 
nominator <- integrate(fc, -Inf, Inf, y=y, eta=eta,delta=delta)$value

# first derivative of the full conditional
f_I <- function(xi, y, eta, delta){
    res <- -delta*(xi-eta) + y - exp(xi)
    return(res)
}

# second derivative of the full conditional
f_II <- function(xi, y, eta, delta){
  res <- -delta - exp(xi)
  return(res)
}

GMRF_approx <- function(xi, xi_0, y, eta, delta){
  a <- fc(xi_0, y, eta, delta) - f_I(xi_0, y, eta, delta)*xi_0 + 0.5*f_II(xi_0, y, eta, delta)*xi_0^2
  b <- f_I(xi_0, y, eta, delta) - f_II(xi_0, y, eta, delta)*xi_0
  c <- -f_II(xi_0, y, eta, delta)
  res <- exp(-1/2 * c* xi^2 + b*xi)
  #norm <- sqrt(delta)/(factorial(y) * sqrt(2*pi))
  norm <- 1
  return(list(res=res*norm, a=a, b=b, c=c))
}

Taylor_approx <- function(xi, xi_0, y, eta, delta){
  res <- fc(xi_0, y, eta, delta) + f_I(xi_0, y, eta, delta)*(xi-xi_0) + 0.5*f_II(xi_0, y, eta, delta)*(xi-xi_0)^2
  #norm <- sqrt(delta)/(factorial(y) * sqrt(2*pi))
  norm <- 1
  return(exp(res)*norm)
}



########################################################################
### EXAMPLE 6.6 HIGH DIMENSIONAL DISTRIBUTION
#########################################################################
# p = number of dimensions
# n = number of points sampled
# N.eff = effective sample size
# X = sample
# w = importance weights
# rejuvs = number of rejuvenations
# normal.draw = sample from the standard normal envelope.
# ftxt = f_t(x_1:t) in book notation
# ft1xt1= f_t-1(x_1:t-1)
# Note: f_t(x_t|x_1:t-1) = ftxt/ft1xt1
#########################################################################

##INITIAL VALUES
p = 50
n = 5000

# sample X_t at time t corresponds to column in matrix X
# the i'th row is X_{1:p}^(i)
X = matrix(NA, nrow = n, ncol = p) 

w = rep(1, n) / n # initialize weights with uniform 1/n
N.eff = rep(NA, p + 1)
N.eff[1] = n # corresponds to uniform weight distribution
ft1xt1 = w
rejuvs = 0

# store weights
ws = matrix(NA, nrow = n, ncol = p)
# target density function
dens.fun = function(x) {
  temp = function(x) {
    exp(-(abs(sqrt(sum(x ** 2))) ** 3) / 3)
  }
  # apply the above density function temp to each row (dimension "1") of matrix x
  # i.e. compute value of density function across all samples
  c(apply(x, 1, temp))
}

# 1D variant of above
dens.fun.1 = function(x) {
  exp(-(abs(sqrt(x ** 2)) ** 3) / 3)
}
ncoeff = integrate(dens.fun.1, -999, 999)$value
dens.fun.1.n = function(x) {
  dens.fun.1(x) / ncoeff
}
##MAIN LOOP
set.seed(4567)
# loop over dimensions / time steps
for (i in 1:p) { 
  # Generate new sample: n points x_t^(i) from the envelope function N(0,1)
  normal.draw = rnorm(n, mean = 0, sd = 1) 
  
  # Compute values of the target density at x_{1:t}
  if (i == 1) {
    # At first time step, compute a one-dimensional distribution
    ftxt = dens.fun(cbind(normal.draw))
  } else {
    # At later time steps, append new sample to previous to make X_t
    Xtemp = cbind(X[, 1:(i - 1)], normal.draw) 
    ftxt = dens.fun(Xtemp)
  }
  # Store new sample to i'th col of matrix X
  X[, i] = normal.draw 
  
  # See comment (*)
  w = w * ftxt / (ft1xt1 * dnorm(normal.draw))
  
  # Store values at current time step to ft1xt1
  ft1xt1 = ftxt
  
  # Normalize and store weights
  w = w / sum(w)
  ws[, i] = w
  
  # Compute effective sample size
  N.eff[i + 1] = 1 / sum(w ** 2)
  
  # If ESS is less than predetermined threshold - here 1/5 of sample size,
  # perform resampling/rejuvenation
  if (N.eff[i + 1] < N.eff[1] / 5) {
    rejuvs = rejuvs + 1
    
    # Sample indices j from (1,...,n) with replacement, with probabilities w_j
    idx = sample(1:n, n, replace = T, prob = w)
    
    # Keep rows corresponding to indices idx from sample matrix X
    X = X[idx,]
    
    # Re-compute density function
    ft1xt1 = dens.fun(X[, 1:i])
    
    # Reset sample weights
    w = rep(1, n) / n
  }
}
plot(0:p, N.eff, type = 'l')

x.avg <- apply(X, 2, mean)
x.ci <- apply(X, 2, quantile, c(0.025, 0.975))
plot(1:p, x.avg, ylim = c(-2, 2))
lines(1:p, x.ci[1, ])
lines(1:p, x.ci[2, ])

x.dens <- density(X)
plot(x.dens, col = 'coral3', lwd = 1.5)
v <- seq(-4, 4, 0.01)
lines(v, dens.fun.1.n(v), col = 'goldenrod3', lwd = 1.5)
lines(v, dnorm(v), col = 'forestgreen', lwd = 1.5)

plot(ws[42, ], type = 'l')

ws.q <- apply(ws, 2, quantile, c(0.025, 0.5, 0.975))
plot(ws.q[2, ], type = 'l', ylim = c(0, 0.0018), col = 'khaki3', lwd = 2)
lines(ws.q[1, ], col = 'hotpink', lwd = 2)
lines(ws.q[3, ], col = 'mistyrose4', lwd = 2)
abline(h = 0.0002, col = 'sandybrown')

plot(density(ws), col = 'tomato')

#NOW TRY TO EXPEND ALL EFFORT ON ONE SIR OF 5000 PTS
# set.seed(11)
x = matrix(rnorm(50 * 5000), 5000, 50)
g = dnorm(x)
f = dens.fun(x)

# Importance weights: Denominator is N_p(0,I_p)
impwt = f / apply(g, 1, prod)
impwt = impwt / sum(impwt)
#ESS
1 / sum(impwt ** 2)  #Answer varies with random number seed

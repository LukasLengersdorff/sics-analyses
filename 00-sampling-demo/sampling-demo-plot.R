############################################################
### HMC sampling demo (Figure 2) ###########################
############################################################

library(tidyverse)
library(sics)
library(rstan)

## compile demo model
mod = stan_model("sampling-demo.stan")

## Parameters for the demo
# Distribution parameters
n = 2                   # dimension (only 2 works for the plots here, but 3 could be interesting to explore with plotly)
mu = c(3,2)             # mean of MV normal
Sigma = diag(c(1,1))    # VCV of MV normal
constraint = LinearConstraint("x1 > x2 > 0.5*x1") # Constraint

# HMC parameters
N = 1000   # Number of samples
dt = 0.01 # HMC time step
nS = 100  # number of HMC leapfrog steps
M = diag(rep(1, n)) # HMC mass matrix
startpoint = mu     # Start point of HMC

seed = 12345

## Prepare SICSFrame-like object with necessary input
SF = list(n = n,
          mu= mu,
          Sigma = Sigma)

## Constraints

cons = list(unc = EmptyConstraint(2), # For unconstrained sampling, panel A
            con = constraint)         # For constrained sampling, panels A and B
res = list()

# Get HMC samples and trajectories from unconstrained and constrained model

for (con_name in names(cons)) {

  set.seed(seed)

  constraint = cons[[con_name]]

  ## Create SICS-like input
  standata = c(SF,
               list(SICS_mA = constraint$mA,
                    SICS_mB = constraint$mB,
                    SICS_mC = constraint$mC,
                    SICS_mU = constraint$mU,
                    SICS_MA = constraint$MA,
                    SICS_MB = constraint$MB,
                    SICS_MC = constraint$MC,
                    SICS_MU = constraint$MU,
                    SICS_a  = as.array(constraint$a),
                    SICS_b  = as.array(constraint$b),
                    SICS_c0 = as.array(constraint$c0),
                    SICS_c1 = as.array(constraint$c1)))

  # Dummy fit to expose STAN's gradient
  fit = suppressWarnings( # To suppress the nearly inevitable 1 divergent transition after "warmup". BTW, I find it weird that rstan does not provide functionality to expose STAN's gradient without the workaround of a one-iteration dummy fit...
    sampling(mod, data = standata, chains = 1, iter = 1, refresh = 0)
    )

  ## --- HMC ---
  leapfrogger = function(x, p, dt, gradient) { # Leapfrog integrator
    p_half = p - 0.5*dt*gradient(x)
    xn = x + dt*p_half
    pn = p_half - 0.5*dt*gradient(xn)
    list(xn = xn, pn = pn)
  }

  grad = function(x) -grad_log_prob(fit, x)    # Gradient
  Minv = solve(M)   # inverse HMC mass matrix

  x0 = transform_par(constraint, startpoint, unbounded= TRUE) # The first point needs to be in the unbounded version of constrained space

  eta = matrix(NA, nrow = N+1, ncol = n)
  eta[1,] = x0

  eta_steps = vector("list", N+1)

  n_same = rep(1, N+1)

  for (i in 2:(N+1)) {
    td = TRUE
    while (td) {
      eta_steps[[i]] = matrix(NA, nrow = nS+1, ncol = n)
      eta_steps[[i]][1,] = eta[i-1,]
      p0 <- p <- mvtnorm::rmvnorm(1, sigma = M)
      for (j in 2:(nS+1)) {
        out = leapfrogger(eta_steps[[i]][j-1,], p, dt, grad)
        eta_steps[[i]][j,] = out$xn
        p = out$pn
      }
      a0 = -log_prob(fit, eta_steps[[i]][1,]) + 0.5*p0 %*% Minv %*% t(p0)
      an = -log_prob(fit, eta_steps[[i]][nS+1,]) + 0.5*p %*% Minv %*% t(p)
      #print(c(an,a0))
      if (runif(1) < exp(an - a0)) {
        eta[i,] = eta_steps[[i]][nS+1,]
        td = FALSE
      } else {
        n_same[i-1] = n_same [i-1]+1
      }
    }
  }

  etas = do.call(rbind, eta_steps)

  res[[con_name]] = list(eta = eta,
                         etas = etas)
}


## Get simulated data in order
res$unc$theta = res$unc$eta
res$unc$thetas = res$unc$etas
res$con$theta = apply(res$con$eta, 1, function(x) transform_par_inv(constraint, x, unbounded = TRUE)) %>% t
res$con$thetas = apply(res$con$etas, 1, function(x) transform_par_inv(constraint, x, unbounded = TRUE)) %>% t
res$con$eta_bounded = apply(res$con$theta, 1, function(x) transform_par(constraint, x, unbounded = FALSE)) %>% t
res$con$etas_bounded = apply(res$con$thetas, 1, function(x) transform_par(constraint, x, unbounded = FALSE)) %>% t

df_u = data.frame(theta1 = res$unc$theta[,1], theta2 = res$unc$theta[,2])
df_us = data.frame(theta1 = res$unc$thetas[,1], theta2 = res$unc$thetas[,2])
df_c = data.frame(theta1 = res$con$theta[,1], theta2 = res$con$theta[,2], eta1 = res$con$eta_bounded[,1], eta2 = res$con$eta_bounded[,2] )
df_cs = data.frame(theta1 = res$con$thetas[,1], theta2 = res$con$thetas[,2], eta1 = res$con$etas_bounded[,1], eta2 = res$con$etas_bounded[,2])
df_us$sample = rep(1:N, each = nS+1)
df_cs$sample = rep(1:N, each = nS+1)

## Density contours
if (TRUE) {
  df_dens = expand.grid(x = seq(0,6,length.out = 400), y = seq(-1, 5, length.out=400))
  f = function(x,y) mvtnorm::dmvnorm(c(x,y), mean = mu, sigma = Sigma)
  f = Vectorize(f, c("x","y"))
  df_dens$dens = with(df_dens, f(x,y))
  df_dens$dens = df_dens$dens/(sum(df_dens$dens))

  # approximate contours
  tmp = data.frame(dens = sort(df_dens$dens))
  tmp$p = cumsum(tmp$dens)
  probs = c(0, 0.25, 0.5, 0.75, 0.95)
  levels = sapply(probs, function(p) min(tmp$dens[tmp$p >= (1-p)]))
  contourfun = Vectorize(function(x) min(levels[levels > x-0.0000001]))
  df_dens$dens_contour = with(df_dens, contourfun(dens))
}

if (TRUE) {
  suppressWarnings({w = sapply(1:nrow(df_dens), function(i) transform_par(constraint, c(df_dens$x[i], df_dens$y[i]), unbounded = FALSE)) %>% t})
  df_dens$w1 = w[,1]
  df_dens$w2 = w[,2]
}

## Plots

orig_limits = list(x = c(0,6),
                   y = c(-1,5))

con_limits = list(x = c(0,2.5),
                  y = c(0,2.5))

psize = 1
palpha = 0.8
pshape = 21
pstroke = 0

lwidth = 0.1
lalpha = 0.5

title_size = 10
title_hjust = 0.5

N_draw = 100

p1 = ggplot() +
  geom_tile(data = df_dens, aes(x = x, y = y, fill = dens_contour)) +
  scale_fill_gradient(low = "white", high = "slateblue", guide = NULL) +
  geom_path(data = df_us[1:(N_draw*(nS+1)),], mapping = aes(x = theta1, y = theta2, group = sample), linewidth = lwidth, alpha = lalpha) +
  geom_point(data = df_u[1:N_draw,], mapping = aes(x = theta1, y = theta2), size = psize, alpha = palpha, shape = pshape, stroke= pstroke, fill = "black") +
  geom_polygon(mapping = aes(x = c(0,6,0,0), y = c(0,6,6,0)), alpha = 0.5, fill = "white") +
  geom_polygon(mapping = aes(x = c(0,6,6,0), y = c(0,3,-1,-1)), alpha = 0.5, fill = "white") +
  geom_abline(slope = 1, intercept = 0) +
  geom_abline(slope = 0.5, intercept = 0) +
  theme_classic() +
  coord_cartesian(xlim = orig_limits$x, ylim = orig_limits$y, expand = FALSE, clip = "on") +
  scale_x_continuous(
    name = expression(theta~""[1]),
    breaks = floor(orig_limits$x[1]):ceiling(orig_limits$x[2])) +
  scale_y_continuous(
    name = expression(theta~""[2]),
    breaks = floor(orig_limits$y[1]):ceiling(orig_limits$y[2])) +
  ggtitle("Unconstrained sampling") +
  theme(plot.title = element_text(size = title_size,  hjust = title_hjust))
#p1

p2 = ggplot() +
  geom_tile(data = df_dens, aes(x = x, y = y, fill = dens_contour)) +
  scale_fill_gradient(low = "white", high = "slateblue", guide = NULL) +
  geom_path(data = df_cs[1:(N_draw*(nS+1)),], mapping = aes(x = theta1, y = theta2, group = sample), linewidth = lwidth, alpha = lalpha) +
  geom_point(data = df_c[1:N_draw,], mapping = aes(x = theta1, y = theta2), size = psize, alpha = palpha, shape = pshape, stroke= pstroke, fill = "black") +
  geom_polygon(mapping = aes(x = c(0,6,0,0), y = c(0,6,6,0)), alpha = 0.5, fill = "white") +
  geom_polygon(mapping = aes(x = c(0,6,6,0), y = c(0,3,-1,-1)), alpha = 0.5, fill = "white") +
  geom_abline(slope = 1, intercept = 0) +
  geom_abline(slope = 0.5, intercept = 0) +
  theme_classic() +
  coord_cartesian(xlim = orig_limits$x, ylim = orig_limits$y, expand = FALSE, clip = "on") +
  scale_x_continuous(
    name = expression(theta~""[1]),
    breaks = floor(orig_limits$x[1]):ceiling(orig_limits$x[2])) +
  scale_y_continuous(
    name = expression(theta~""[2]),
    breaks = floor(orig_limits$y[1]):ceiling(orig_limits$y[2])) +
  ggtitle("SICS: original space") +
  theme(plot.title = element_text(size = title_size,  hjust = title_hjust))
#p2

p3 = ggplot() +
  geom_tile(data = df_dens, aes(x = w1, y = w2, fill = dens_contour)) +
  scale_fill_gradient(low = "white", high = "slateblue", guide = NULL) +
  geom_path(data = df_cs[1:(N_draw*(nS+1)),], mapping = aes(x = eta1, y = eta2, group = sample), linewidth = lwidth, alpha = lalpha) +
  geom_point(data = df_c[1:N_draw,], mapping = aes(x = eta1, y = eta2), size = psize, alpha = palpha, shape = pshape, stroke= pstroke, fill = "black") +
  theme_classic() +
  coord_cartesian(xlim = con_limits$x, ylim = con_limits$y, expand = FALSE, clip = "on") +
  scale_x_continuous(
    name = expression(omega~~""[1]),
    breaks = floor(con_limits$x[1]):ceiling(con_limits$x[2])) +
  scale_y_continuous(
    name = expression(omega~~""[2]),
    breaks = floor(con_limits$x[1]):ceiling(con_limits$x[2])) +
  ggtitle("SICS: transformed space") +
  theme(plot.title = element_text(size = title_size,  hjust = title_hjust))
#p3

P = ggpubr::ggarrange(p1,p2,p3, nrow = 1, labels = c("A","B","C"))
P

ggsave("sampling-demo-plot.png", P, width = 21, height = 21/3 + 0.5, dpi = 300, units = "cm")

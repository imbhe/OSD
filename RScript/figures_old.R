
# Figure 1.
mu <- opt_sampling_scheme(n, G, H, w, crit = "L", L = sqrtm(H), "PO-WR")$mu

fig1 <-  ggplot(data, aes(x = impact_speed0)) + 
  stat_function(aes(lty = "Full data", colour = "Full data"), fun = dlnorm, args = list(meanlog = theta0[1], sdlog = theta0[2])) +
  labs(x = "Baseline impact speed (km/h)", 
       y = "Density", 
       lty = NULL, 
       colour = NULL) + 
  theme(legend.position = c(0.95, 0.95), 
        legend.justification = c("right", "top"))

for ( i in 1:10) {
  wt <- w * rpois(n = N, lambda = mu) / mu
  ix <- which(wt > 0)
  thetahat <- c(weighted.mean(logy[ix], wt[ix]), sqrt(cov.wt(matrix(logy[ix], ncol = 1), wt[ix], method = "ML")$cov)) 
  fig1 <- fig1 + stat_function(aes(lty = "Subsample", colour = "Subsample"), fun = dlnorm, args = list(meanlog = thetahat[1], sdlog = thetahat[2])) 
}

fig1 <- fig1 + 
  stat_function(aes(lty = "Full data", colour = "Full data"), fun = dlnorm, args = list(meanlog = theta0[1], sdlog = theta0[2])) +
  scale_colour_manual(values = c("black", "grey70"))
ggsave("Figures/Figure1.png", fig1, dpi = 1000, width = 80, height = 50, unit = "mm")

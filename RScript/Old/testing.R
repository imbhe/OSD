# Functions.
theta <- function(Y, w = rep(1, nrow(Y))) {
  colSums(w * Y) / sum(w)
}

grads <- function(theta0, Y, w = rep(1, nrow(Y))) {
  return(t(-w * t(t(Y) - theta0)))
}

hess <- function(theta0, w = rep(1, nrow(Y))) {
  return(diag(sum(w), length(theta0)))
}

# Data.
data <- VSdata %>% 
  filter(crash0 == 1) 

N <- nrow(data)
n <- max(round(N / 100), 10)
w <- with(data, crash0 * prob)

# As matrix.
Y <- data %>% 
  dplyr::select(impact_speed_reduction, injury_risk_reduction, crash_avoidance) %>%
  as.matrix()
(theta0 <- theta(Y, w))
G <- grads(theta0, Y, w)
H <- hess(theta0, w)
# Scaled.
sds <- sqrt(diag(cov.wt(Y, w, method = "ML")$cov))
A <- diag(1 / sds)
Y_scal <- Y %*% t(A)
theta0_scal <- theta(Y_scal, w)
G_scal <- grads(theta0_scal, Y_scal, w)
H_scal <- hess(theta0_scal, w)
cov.wt(Y_scal, w, method = "ML")$cov

# Centered.
Y_cent <- t(t(Y) - theta0)
theta0_cent <- theta(Y_cent, w)
G_cent <- grads(theta0_cent, Y_cent, w)
H_cent <- hess(theta0_cent, w)

# Scaled and centered.
sds <- sqrt(diag(cov.wt(Y, w, method = "ML")$cov))
A <- diag(1 / sds)
Y_scal_cent <- t(t(Y) - theta0) %*% t(A)
theta0_scal_cent <- theta(Y_scal_cent, w)
G_scal_cent <- grads(theta0_scal_cent, Y_scal_cent, w)
H_scal_cent <- hess(theta0_scal_cent, w)

# Orthogonalised.
B <- solve(sqrtm(cov.wt(Y, w, method = "ML")$cov))
Y_orth <- Y %*% t(B)
theta0_orth <- theta(Y_orth, w)
G_orth <- grads(theta0_orth, Y_orth, w)
H_orth <- hess(theta0_orth, w)

# Arbitrary affine transformation.
a <- rnorm(3)
C <- matrix(runif(9), 3)
C <- C / rowSums(C)
Y_trans <- t(a + t(Y_scal_cent %*% C))
theta0_trans <- theta(Y_trans, w)
G_trans <- grads(theta0_trans, Y_trans, w)
H_trans <- hess(theta0_trans, w)


# D-optimality on original and transformed data.
res1 <- opt_sampling_scheme(n, G,           H,           w, crit = "D")
res2 <- opt_sampling_scheme(n, G_cent,      H_cent,      w, crit = "D")
res3 <- opt_sampling_scheme(n, G_scal,      H_scal,      w, crit = "D")
res4 <- opt_sampling_scheme(n, G_scal_cent, H_scal_cent, w, crit = "D")
res5 <- opt_sampling_scheme(n, G_orth,      H_orth,      w, crit = "D")
res6 <- opt_sampling_scheme(n, G_trans,     H_trans,     w, crit = "D")
par(mfrow = c(2, 3))
plot(res1$mu, res2$mu, bty = "l")
plot(res1$mu, res3$mu, bty = "l")
plot(res1$mu, res4$mu, bty = "l")
plot(res1$mu, res5$mu, bty = "l")
plot(res1$mu, res6$mu, bty = "l")
par(mfrow = c(1, 1))

res7 <- opt_sampling_scheme(n, G_cent, H_cent,  w, crit = "L", L = sqrtm(solve(acov(res1$mu, G, H, w))), init = res1$mu)


# Now try PO-WOR.
res11 <- opt_sampling_scheme(n, G,           H,           w, crit = "D", design = "PO-WOR")
res12 <- opt_sampling_scheme(n, G_cent,      H_cent,      w, crit = "D", design = "PO-WOR")
res13 <- opt_sampling_scheme(n, G_scal,      H_scal,      w, crit = "D", design = "PO-WOR")
res14 <- opt_sampling_scheme(n, G_scal_cent, H_scal_cent, w, crit = "D", design = "PO-WOR")
res15 <- opt_sampling_scheme(n, G_orth,      H_orth,      w, crit = "D", design = "PO-WOR")
res16 <- opt_sampling_scheme(n, G_trans,     H_trans,     w, crit = "D", design = "PO-WOR")
par(mfrow = c(2, 3))
plot(res11$mu, res12$mu, bty = "l")
plot(res11$mu, res13$mu, bty = "l")
plot(res11$mu, res14$mu, bty = "l")
plot(res11$mu, res15$mu, bty = "l")
plot(res11$mu, res16$mu, bty = "l")
par(mfrow = c(1, 1))

plot(res1$mu, res11$mu, bty = "l")
res11 <- opt_sampling_scheme(n, G, H, w, crit = "D", design = "PO-WOR", init = res1$mu)


# E-optimality on original and transformed data.
# Invariant to coordinate-wise location shifts. Not invariant otherwise.
# Converges for non-standardised data because there is one clear dominant direction. 
# Does not converge otherwise.
res1 <- opt_sampling_scheme(n, G,           H,           w, crit = "E")
res2 <- opt_sampling_scheme(n, G_cent,      H_cent,      w, crit = "E")
res3 <- opt_sampling_scheme(n, G_scal,      H_scal,      w, crit = "E")
res4 <- opt_sampling_scheme(n, G_scal_cent, H_scal_cent, w, crit = "E")
res5 <- opt_sampling_scheme(n, G_orth,      H_orth,      w, crit = "E")
res6 <- opt_sampling_scheme(n, G_trans,     H_trans,     w, crit = "E")


# Approximate E-optimality using Phi_p optimality.
res1 <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "Phi_p", p = 1)
res2 <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "Phi_p", p = 2)
res3 <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "Phi_p", p = 3)
res4 <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "Phi_p", p = 4)
res5 <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "Phi_p", p = 5, init = res4$mu)
res6 <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "Phi_p", p = 6, init = res5$mu)
res7 <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "Phi_p", p = 7, init = res6$mu)
res8 <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "Phi_p", p = 8, init = res7$mu)
Phi(acov(res1$mu, G_scal, H_scal), "E") / Phi(acov(res7$mu, G_scal, H_scal), "E") # 8% improvement.

res <- approx_E_opt(n, G_scal, H_scal, w)
Phi(acov(res$mu, G_scal, H_scal), "E") / Phi(acov(res1$mu, G_scal, H_scal), "E") # 8% improvement compared to A-optimality.



# A-optimality on original and transformed data.
# Invariant to coordinate-wise location shifts. 
# Not invariant to arbitrary affine transformations.
res1 <- opt_sampling_scheme(n, G,       H,       w, crit = "A")
res2 <- opt_sampling_scheme(n, G_scal,  H_scal,  w, crit = "A")
res3 <- opt_sampling_scheme(n, G_cent,  H_cent,  w, crit = "A")
res4 <- opt_sampling_scheme(n, G_orth,  H_orth,  w, crit = "A")
res5 <- opt_sampling_scheme(n, G_trans, H_trans, w, crit = "A")
res6 <- opt_sampling_scheme(n, G,       H,       w, crit = "Phi_p", r = 1)
par(mfrow = c(2, 2))
plot(res1$mu, res2$mu, bty = "l")
plot(res1$mu, res3$mu, bty = "l")
plot(res1$mu, res4$mu, bty = "l")
plot(res1$mu, res5$mu, bty = "l")
par(mfrow = c(1, 1))

plot(res1$mu, res6$mu, bty = "l")


# L-optimality wrt. matrix square root/Cholesky decomposition of information matrix. 
# In this case equivalent to A-optimality.
L <- t(chol(H))
res <- opt_sampling_scheme(n, G, H, w, crit = "L", L = L, kappa_tol = 1e7)


# L-optimality wrt. matrix square root/Cholesky decomposition of empirical covariance matrix of data. 
# Invariant to coordinate-wise affine transformations.
L <- t(chol(cov.wt(Y, w, method = "ML")$cov))
res1 <- opt_sampling_scheme(n, G,      H,      w, crit = "L", L = L, kappa_tol = 1e7)
res2 <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "L", L = L, kappa_tol = 1e7)
res3 <- opt_sampling_scheme(n, G_cent, H_cent, w, crit = "L", L = L, kappa_tol = 1e7)
res4 <- opt_sampling_scheme(n, G_orth, H_orth, w, crit = "L", L = L, kappa_tol = 1e7)
res5 <- opt_sampling_scheme(n, G_trans, H_trans, w, crit = "L", L = L, kappa_tol = 1e7)
res6 <- opt_sampling_scheme(n, G_trans2, H_trans2, w, crit = "L", L = L, kappa_tol = 1e7)
par(mfrow = c(2, 3))
plot(res1$mu, res2$mu, bty = "l")
plot(res1$mu, res3$mu, bty = "l")
plot(res1$mu, res4$mu, bty = "l")
plot(res1$mu, res5$mu, bty = "l")
plot(res1$mu, res6$mu, bty = "l")
par(mfrow = c(1, 1))


# L-optimality wrt. matrix square root/Cholesky decomposition of full-data empirical covariance matrix. 
# Invariant to arbitrary affine transformations! 
# In this case equivalent to L-optimality wrt. matrix square root/Cholesky decomposition of empirical covariance matrix of data. 
GtG <- matrix(0, nrow(G), nrow(G))
for ( i in 1:ncol(G) ) {
  GtG <- GtG + tcrossprod(G[, i])
}
L <- solve(H, t(chol(GtG)))
res1 <- opt_sampling_scheme(n, G,      H,      w, crit = "L", L = L, kappa_tol = 1e7)
res2 <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "L", L = L, kappa_tol = 1e7)
res3 <- opt_sampling_scheme(n, G_cent, H_cent, w, crit = "L", L = L, kappa_tol = 1e7)
res4 <- opt_sampling_scheme(n, G_orth, H_orth, w, crit = "L", L = L, kappa_tol = 1e7)
res5 <- opt_sampling_scheme(n, G_trans, H_trans, w, crit = "L", L = L, kappa_tol = 1e7)
res6 <- opt_sampling_scheme(n, G_trans2, H_trans2, w, crit = "L", L = L, kappa_tol = 1e7)
par(mfrow = c(2, 3))
plot(res1$mu, res2$mu, bty = "l")
plot(res1$mu, res3$mu, bty = "l")
plot(res1$mu, res4$mu, bty = "l")
plot(res1$mu, res5$mu, bty = "l")
plot(res1$mu, res6$mu, bty = "l")
par(mfrow = c(1, 1))


# c-optimality.
res1 <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "L", L = c(1, 0, 0), design = "PO-WOR", kappa_tol = 1e6)
res2 <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "L", L = c(1, 0, 0), design = "PO-WR" , kappa_tol = 1e6)
res3 <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "L", L = c(0, 1, 0), design = "PO-WOR", kappa_tol = 1e6)
res4 <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "L", L = c(0, 1, 0), design = "PO-WR" , kappa_tol = 1e6)
res5 <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "L", L = c(0, 0, 1), design = "PO-WOR", kappa_tol = 1e6)
res6 <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "L", L = c(0, 0, 1), design = "PO-WR" , kappa_tol = 1e6)

# c-optimality on standardised and orthogonalised data.
# Invariant under coordinate-wise affine transformations.
# Invariant under arbitrary affine transformations if the linear combination of interest is transformed along with the data.
L <- c(0, 1, 0)
res7  <- opt_sampling_scheme(n, G_scal, H_scal, w, crit = "L", L = L, design = "PO-WR")
res8  <- opt_sampling_scheme(n, G_cent, H_cent, w, crit = "L", L = L, design = "PO-WR", kappa_tol = 1e7)
res9  <- opt_sampling_scheme(n, G_orth, H_orth, w, crit = "L", L = L, design = "PO-WR")
res10 <- opt_sampling_scheme(n, G_orth, H_orth, w, crit = "L", L = as.numeric(A %*% L), design = "PO-WR")
par(mfrow = c(2, 2))
plot(res4$mu, res7$mu, bty = "l")
plot(res4$mu, res8$mu, bty = "l")
plot(res4$mu, res9$mu, bty = "l")
plot(res4$mu, res10$mu, bty = "l")
par(mfrow = c(1, 1))



# D- and E-optimality.
res1 <- opt_sampling_scheme(n, G, H, w, crit = "D", design = "PO-WOR")
res2 <- opt_sampling_scheme(n, G, H, w, crit = "D", design = "PO-WR" )
res3 <- opt_sampling_scheme(n, G, H, w, crit = "E", design = "PO-WOR")
res4 <- opt_sampling_scheme(n, G, H, w, crit = "E", design = "PO-WR" )

# D- and E-optimality, start from random scheme.
res1 <- opt_sampling_scheme(n, G, H, w, crit = "D", design = "PO-WOR", init = runif(N))
res2 <- opt_sampling_scheme(n, G, H, w, crit = "D", design = "PO-WR",  init = runif(N))
res3 <- opt_sampling_scheme(n, G, H, w, crit = "E", design = "PO-WOR", init = runif(N))
res4 <- opt_sampling_scheme(n, G, H, w, crit = "E", design = "PO-WR" , init = runif(N))

# D- and E-optimality, start from density sampling scheme. 
res1 <- opt_sampling_scheme(n, G, H, w, crit = "D", design = "PO-WOR", init = data$prob)
res2 <- opt_sampling_scheme(n, G, H, w, crit = "D", design = "PO-WR" , init = data$prob)
res3 <- opt_sampling_scheme(n, G, H, w, crit = "E", design = "PO-WOR", init = data$prob)
res4 <- opt_sampling_scheme(n, G, H, w, crit = "E", design = "PO-WR" , init = data$prob)
plot(res2$mu, res1$mu, bty = "l")


# A-optimality
res1 <- opt_sampling_scheme(n, G, H, w, design = "PO-WOR", crit = "A")
res2 <- opt_sampling_scheme(n, G, H, w, design = "PO-WR",  crit = "A")
plot(res2$mu, res1$mu)


# L-optimality wrt. matrix square root/Cholesky decomposition of information matrix. 
L <- t(chol(H))
res1 <- opt_sampling_scheme(n, G, H, w, crit = "L", L = L, design = "PO-WOR")
res2 <- opt_sampling_scheme(n, G, H, w, crit = "L", L = L, design = "PO-WR")


# L-optimality wrt. matrix square root/Cholesky decomposition of empirical precision matrix. 
GtG <- matrix(0, nrow(G), nrow(G))
for ( i in 1:ncol(G) ) {
  GtG <- GtG + tcrossprod(G[, i])
}
L <- solve(H, t(chol(GtG)))
res3 <- opt_sampling_scheme(n, G, H, w, crit = "L", L = L, design = "PO-WOR")
res4 <- opt_sampling_scheme(n, G, H, w, crit = "L", L = L, design = "PO-WR")

plot(res1$mu, res3$mu)


# c-optimality
res1 <- opt_sampling_scheme(n, G, H, w, crit = "L", L = c(1, 0), design = "PO-WOR")
res2 <- opt_sampling_scheme(n, G, H, w, crit = "L", L = c(1, 0), design = "PO-WR")
res3 <- opt_sampling_scheme(n, G, H, w, crit = "L", L = c(0, 1), design = "PO-WOR")
res4 <- opt_sampling_scheme(n, G, H, w, crit = "L", L = c(0, 1), design = "PO-WR")





# %>% 
#   filter(!(caseID %in% c(23, 25, 41)))

for ( i in 1:length(unique(VSdata$caseID)) ) {
  id <- unique(VSdata$caseID)[i]
  
  dat <- VSdata %>%  
    filter(caseID == id)
  
  y <- with(dat, impact_speed0 / max_impact_speed)
  ymax <- dat$max_impact_speed
  w <- with(dat, prob)
  fit <- glm(y~dec*OEOFF, data = dat, family = "quasibinomial")
  pred <- ymax * fit$fitted
  
  plt <- dat %>% mutate(impact_speed0  = as.numeric(pred), grp = 0) %>%
    bind_rows(dat %>% mutate(grp = 1))
  
  p <- ggplot(plt) +
    geom_line(aes(x = OEOFF, y = impact_speed0, group = dec, colour = dec)) +
    facet_wrap(~grp) + 
    labs(title = id) + 
    theme_classic() +
    theme(legend.position = NULL)
  
  print(p)
  
}

data$caseID <- factor(data$caseID)
data$y <- with(data, impact_speed0 / max_impact_speed)

fit <- glm(y~caseID*dec*OEOFF, data = data, family = "quasibinomial")
pred <- predict(fit, type = "response")

plt <- data %>% 
  mutate(impact_speed  = max_impact_speed * as.numeric(pred), grp = 1) %>%
  bind_rows(data %>% mutate(impact_speed = impact_speed0, grp = 0))


for ( i in 1:1 ){
  
  p <- ggplot(plt %>% filter(caseID %in% (4 * (i - 1) + 1):(4 * i)) ) +
    geom_line(aes(x = OEOFF, y = impact_speed, group = dec, colour = dec))  +
    facet_grid(caseID~grp, labeller = labeller(grp = c("0" = "Observed", "1" = "Predicted"))) +
    theme(legend.position = NULL)
  
  print(p)
}


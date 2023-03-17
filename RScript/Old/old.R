



# Sampling schemes.
params <- tibble(crit = c("A", "D", "E", rep("c", 3), rep("L", 3), rep("Phi_r", 3)), 
                 name = c("A", "D", "E", paste0("c", 1:3), paste0("L", 1:3), "Phi0.5", "Phi5", "Phi10"),
                 p = c(rep(NA, 9), 0.5, 5, 10))
L <- list(NULL, NULL, NULL, c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), sqrtm(H), sqrtm(H %*% solve(V)), sqrtm(solve(cov.wt(Y, w)$cov)), NULL, NULL, NULL)
params$L <- L
params <- crossing(design = c("PO-WR", "PO-WOR"), params, niter = NA, status = NA, A = NA, D = NA, E = NA, c1 = NA, c2 = NA, c3 = NA, L1 = NA, L2 = NA, L3 = NA, Phi05 = NA, Phi5 = NA, Phi10 = NA)

# Find optimal sampling schemes and evaluate performance analytically. 
for ( i in 1:nrow(params) ) { 
  opt <- opt_sampling_scheme(n, G, H, w, crit = params$crit[i], L = params$L[[i]], p = params$p[i], design = params$design[i])
  params$niter[i] <- opt$niter
  params$status[i] <- opt$status
  Gamma <- acov(opt$mu, G, H, w, design = params$design[i])
  params$A[i] <- Phi(Gamma, crit = "A")
  params$D[i] <- Phi(Gamma, crit = "D")
  params$E[i] <- Phi(Gamma, crit = "E")
  params$c1[i] <- Phi(Gamma, crit = "c", L = L[[4]])
  params$c2[i] <- Phi(Gamma, crit = "c", L = L[[5]])
  params$c3[i] <- Phi(Gamma, crit = "c", L = L[[6]])
  params$L1[i] <- Phi(Gamma, crit = "L", L = L[[7]])
  params$L2[i] <- Phi(Gamma, crit = "L", L = L[[8]])
  params$L3[i] <- Phi(Gamma, crit = "L", L = L[[9]])
  params$Phi05[i] <- Phi(Gamma, crit = "Phi_r", r = 0.5)
  params$Phi5[i] <- Phi(Gamma, crit = "Phi_r", r = 5)
  params$Phi10[i] <- Phi(Gamma, crit = "Phi_r", r = 10)
}

# Relative efficiencies.
res <- params %>% 
  group_by(design) %>% 
  mutate(A = min(A) / A, 
         c1 = min(c1) / c1, 
         c2 = min(c2) / c2, 
         c3 = min(c3) / c3,
         D = exp((min(D) - D) / length(theta0)),
         E = min(E) / E,
         L1 = min(L1) / L1,
         L2 = min(L2) / L2,
         L3 = min(L3) / L3,
         Phi05 = min(Phi05) / Phi05,
         Phi5 = min(Phi5) / Phi5,
         Phi10 = min(Phi10) / Phi10) %>% 
  mutate(niter_s = ifelse(status == 0, sprintf("%.0f", niter), "Did not converge")) %>% 
  ungroup() %>%
  dplyr::select(design, crit, name, niter, niter_s, everything())

# If optimal sampling scheme could not be found.
for ( i in 1:nrow(res) ) {
  if ( res$niter_s[i] == "Did not converge") {
    res[i, 9:ncol(res)] <- NA
    if ( res$crit[i] == "E" ) {
      res$E[which(res$design == res$design[i])] <- NA
    } else if ( res$crit[i] == "Phi_r" & res$p[i] == 5 ) {
      res$Phi5[which(res$design == res$design[i])] <- NA
    } else if ( res$crit[i] == "Phi_r" & res$p[i] == 10 ) {
      res$Phi10[which(res$design == res$design[i])] <- NA
    }
  }
}

save(res, file = "Results/tab4.RData")


# Write xtable to file.
load(file = "Results/tab4.RData")
res %>% 
  filter(design == "PO-WR") %>% 
  dplyr::select(name, niter_s, A, c1, c2, c3, D, E) %>% 
  dplyr::rename("Optimality criterion" = name,
                "Number of iterations" = niter_s) %>% 
  xtable(digits = 2, caption = "", label = "tab:finite_population_inference", align = c("llrrrrrrr")) %>% 
  print(caption.placement = "top", 
        table.placement = "ht!", 
        include.rownames = FALSE, 
        NA.string="NA",
        file = "Output/tab4.txt")# Hessian.
# p[p > (1-1e-12)] <- 1-1e-12
# H <- matrix(0, length(theta0), length(theta0))
# for ( i in 1:nrow(X) ) {
#   H <- H + p[i] * (1 - p[i]) * tcrossprod(X[i, ])
# }

X <- model.matrix(~ caseID + dec * OEOFF + dec * I(OEOFF^4), data = data)
ymax <- dat$max_impact_speed
y <- data$impact_speed0

f <- function(beta) {
  yhat <- ymax * (1 + exp(- X %*% beta))^(-1)
  SS <- with(dat, sum((impact_speed - yhat)^2))
  return(SS)
}

fit2 <- optim(rep(0, ncol(X)), f, method = "BFGS", control = list(trace = TRUE))
pred <- ymax * (1 + exp(- X %*% fit$par))^(-1)


V <- matrix(0, 3, 3)
for ( i in 1:(N - 1) ) {
  V <- V + tcrossprod(G[, i])    
}

G[, N] <- rep(0.01, 3)
eps <- c(10^(-(0:15)))
x <- rep(NA, length = length(eps))
for ( i in seq_along(eps) ) {
  mu <- c(rep( (n - eps) / (N - 1), N - 1), eps)
  V_eps <- V / mu[1] + 1 / eps[i] * tcrossprod(G[, N])
  x[i] <- t(G[, N]) %*% solve(V_eps) %*% G[, N] / eps[i]
}
plot(eps, x)
max(x)

# OLS formulation.
data_long <- data %>% 
  gather(impact_speed_reduction, injury_risk_reduction, crash_avoidance, key = "variable", value = y) %>% 
  arrange(caseID, desc(dec), OEOFF) %>% 
  mutate(w = crash0 * prob) %>% 
  dplyr::select(caseID, OEOFF, dec, crash0, prob, w, variable, y)
y <- data_long$y

tic()
mod <- lm(y~-1+variable, data = data_long, weights = w)
toc()
theta <- coef(mod)[c("variableimpact_speed_reduction", 
                     "variableinjury_risk_reduction",
                     "variablecrash_avoidance")]

names(theta) <- c("mean_impact_speed_reduction", 
                  "mean_injury_risk_reduction", 
                  "crash_avoidance_rate")

# Target parameter, basic formulation.
tic()
mean_impact_speed_reduction <- with(data, sum(crash0 *prob * impact_speed_reduction) / sum(crash0 * prob))
mean_injury_risk_reduction <- with(data, sum(crash0 * prob * injury_risk_reduction) / sum(crash0 * prob))
crash_avoidance_rate <- with(data, sum(crash0 * prob * crash_avoidance) / sum(crash0 * prob))
toc()



# constrOptim. ----
# Does not work well in high dimensions.

U1 <- Diagonal(N-1, 1)
c1 <- rep(0, N-1)
U2 <- Diagonal(N-1, -1)
c2 <- rep(-1, N-1) # -1 for PO-WOR, -n otherwise.
U3 <- rep(-1, N-1)
c3 <- -n
ui <- rbind(U1, U2, U3)
ci <- c(c1, c2, c3)

tic()
x0 <- matrix(rep(n / N, N)[-1], nrow = N-1)
opt2 <- constrOptim(x0, 
                    f = function(x) { x <- c(n - sum(x), x) 
                    sum(c^2 / x)
                    }, 
                    grad = NULL,
                    method = "Nelder-Mead",
                    ui, ci)
opt2$par <- c(n - sum(opt2$par), opt2$par)
toc()

set.seed(123)
myopt1 <- optDesign(inclusionprobabilities(runif(N), n), theta, Y, w, "PO-WR", crit = "D") 
plot(myopt1$mu, myopt1$ci)
abline(0, mean(myopt1$ci / myopt1$mu))

myopt2 <- optDesign(inclusionprobabilities(runif(N), n), theta, Y, w, "PO-WR", crit = "D") 
plot(myopt2$mu, myopt2$ci)
abline(0, mean(myopt2$ci / myopt2$mu))

myopt3 <- optDesign(inclusionprobabilities(runif(N), n), theta, Y, w, "PO-WR", crit = "D") 
plot(myopt3$mu, myopt3$ci)
abline(0, mean(myopt3$ci / myopt3$mu))

myopt4 <- optDesign(inclusionprobabilities(runif(N), n), theta, Y, w, "PO-WR", crit = "D") 
plot(myopt4$mu, myopt4$ci)
abline(0, mean(myopt4$ci / myopt4$mu))

myopt5 <- optDesign(inclusionprobabilities(runif(N), n), theta, Y, w, "PO-WR", crit = "D", maxiter = 10) 
plot(myopt5$mu, myopt5$ci)
abline(0, mean(myopt5$ci / myopt5$mu))

myopt6 <- optDesign(inclusionprobabilities(runif(N), n), theta, Y, w, "PO-WR", crit = "D") 
plot(myopt6$mu, myopt6$ci)
abline(0, mean(myopt6$ci / myopt6$mu))

plot(myopt1$mu, myopt2$mu)
abline(0, 1)
plot(myopt1$mu, myopt3$mu)
abline(0, 1)
plot(myopt1$mu, myopt4$mu)
abline(0, 1)
plot(myopt1$mu, myopt5$mu)
abline(0, 1)
plot(myopt1$mu, myopt6$mu)
abline(0, 1)
# All find same solution.



# NLoptr options.
opts <- list("algorithm" = "NLOPT_LD_AUGLAG_EQ",
             "print_level" = 0, 
             "xtol_rel" = 1e-4,
             "local_opts" = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1e-4))

opts <- list("algorithm" = "NLOPT_LN_AUGLAG_EQ",
             "print_level" = 0,
             "xtol_rel" = 1e-4,
             "local_opts" = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-4))
opt <- nloptr( x0 = mu0,
               eval_f = function(x) Phi(vcov(x, theta, Y, w, design = "PO-WOR"), crit = "D"),
               # eval_grad_f = function(x) dPhidMu(x, theta, Y, w, "PO-WOR", crit = "D"),
               # eval_g_eq = function(x) sum(x) - n, 
               # eval_jac_g_eq = function(x) rep(1, length(x)), 
               lb = rep(1e-16, N), 
               ub = rep(1, N), # 1 for PO-WOR, n otherwise.
               opts = opts)

Phi(vcov(mu0, theta, Y, w, design = "PO-WOR"), crit = "D")
opt$objective


psopt <- pso::psoptim(par = mu0, 
                      fn = function(x) Phi(x * n / sum(x), theta, Y, w, "PO-WR", crit = "D"), 
                      lower = 1e-12, 
                      upper = 1, 
                      control = list(trace = 1, 
                                     maxit = 30,
                                     s = 500))
psopt$par <- n * psopt$par / sum(psopt$par)
Phi(psopt$par, theta, Y, w, "PO-WR", crit = "D")
plot(mu0, psopt$par)


opt <- nloptr( x0 = mu1,
               eval_f = function(x) Phi(vcov(x, theta, Y, w, design = "PO-WOR"), crit = "D"),
               # eval_grad_f = function(x) dPhidMu(x, theta, Y, w, "PO-WOR", crit = "D"),
               # eval_g_eq = function(x) sum(x) - n, 
               # eval_jac_g_eq = function(x) rep(1, length(x)), 
               lb = rep(1e-16, N), 
               ub = rep(1, N), # 1 for PO-WOR, n otherwise.
               opts = opts)

Phi(vcov(mu1, theta, Y, w, design = "PO-WOR"), crit = "D")
opt$objective



psopt <- pso::psoptim(par = mu1, 
                      fn = function(x) Phi(x * n / sum(x), theta, Y, w, "PO-WR", crit = "D"), 
                      lower = 1e-12, 
                      upper = 1, 
                      control = list(trace = 1, 
                                     maxit = 30,
                                     s = 500))
psopt$par <- n * psopt$par / sum(psopt$par)
Phi(psopt$par, theta, Y, w, "PO-WR", crit = "D")
plot(mu1, psopt$par)

opt <- nloptr( x0 = mu2,
               eval_f = function(x) Phi(vcov(x, theta, Y, w, design = "PO-WOR"), crit = "D"),
               # eval_grad_f = function(x) dPhidMu(x, theta, Y, w, "PO-WOR", crit = "D"),
               # eval_g_eq = function(x) sum(x) - n, 
               # eval_jac_g_eq = function(x) rep(1, length(x)), 
               lb = rep(1e-16, N), 
               ub = rep(1, N), # 1 for PO-WOR, n otherwise.
               opts = opts)

Phi(vcov(mu2, theta, Y, w, design = "PO-WOR"), crit = "D")
opt$objective


100 * (log(det(Gamma(opt$solution, theta, Y, w, "PO-WOR"))) / log(det(Gamma(mu2, theta, Y, w, "PO-WOR"))) - 1)
log(det(Gamma(opt$solution, theta, Y, w, "PO-WOR")))

psopt <- pso::psoptim(par = mu2, 
                      fn = function(x) Phi(x * n / sum(x), theta, Y, w, "PO-WR", crit = "D"), 
                      lower = 1e-12, 
                      upper = 1, 
                      control = list(trace = 1, 
                                     maxit = 30,
                                     s = 500))
psopt$par <- n * psopt$par / sum(psopt$par)
Phi(mu2, theta, Y, w, "PO-WR", crit = "D")
Phi(psopt$par, theta, Y, w, "PO-WR", crit = "D")
plot(mu2, psopt$par)



opt <- nloptr( x0 = mu4,
               eval_f = function(x) Phi(vcov(x, theta, Y, w, design = "PO-WOR"), crit = "D"),
               # eval_grad_f = function(x) dPhidMu(x, theta, Y, w, "PO-WOR", crit = "D"),
               # eval_g_eq = function(x) sum(x) - n, 
               # eval_jac_g_eq = function(x) rep(1, length(x)), 
               lb = rep(1e-16, N), 
               ub = rep(1, N), # 1 for PO-WOR, n otherwise.
               opts = opts)

Phi(vcov(mu3, theta, Y, w, design = "PO-WOR"), crit = "D")
opt$objective

myopt <- optDesign(mu3, theta, Y, w, design = "PO-WOR", crit = "D")
myopt$val
plot(myopt$mu, myopt$ci)
abline(0, myopt$ci[1] / myopt$mu)


psopt <- pso::psoptim(par = mu3, 
                      fn = function(x) Phi(x * n / sum(x), theta, Y, w, "PO-WR", crit = "D"), 
                      lower = 1e-12, 
                      upper = 1, 
                      control = list(trace = 1, 
                                     maxit = 30,
                                     s = 500))
psopt$par <- n * psopt$par / sum(psopt$par)
Phi(mu3, theta, Y, w, "PO-WR", crit = "D")
Phi(psopt$par, theta, Y, w, "PO-WR", crit = "D")
plot(mu3, psopt$par)



# Derivative of of objective function.
# grad <- function(L, Hinv, mu, theta, Y, w) {
#   
return(vapply(1:N, function(ix) crossprod(w[ix] * L %*% Hinv %*% (theta - Y[ix, ]) / mu[ix]), numeric(1)))
# 
# }



opt1 <- optim(rep(0, 3), function(theta) f(theta, Y, w), method = "BFGS", hessian = TRUE)
opt2 <- optim(rep(0, 3), function(theta) f(theta, Y, w), gr = function(theta) rowSums(grads(theta, Y, w)), method = "BFGS", hessian = TRUE)
hess(theta, Y, w)
opt2$hessian



# Empirical risk function.
f <- function(theta, Y, w) { 
  if ( length(w) != nrow(Y) ) { stop("Dimensions do not agree.")}
  if ( length(theta) != ncol(Y) ) { stop("Dimensions do not agree.")}
  
  N <- nrow(Y)
  losses <- vapply(1:N, function(ix) 0.5 * w[ix] * crossprod(theta - Y[ix, ]), numeric(1))
  
  return(sum(losses))
}



grad <- function(grads) {
  return(rowSums(grads))
}



acov(rep(n / N, N), grads(theta, Y, w), hess(theta, w), w, design = "PO-WR")
acov(rep(n / N, N), grads(theta, Y, w), hess(theta, w), w, design = "PO-WOR")
acov(rep(n / N, N), grads(theta, Y, w), hess(theta, w), w, design = "MULT")


opt6a <- calculate_sampling_scheme(rep(n / N, N), G, H, w, design, crit = "L", L = L)
cov <- cov.wt(Y, wt = w, method = "ML")$cov
L <- t(chol(cov))
opt6b <- calculate_sampling_scheme(rep(n / N, N), G, H, w, design, crit = "L", L = L)


f <- function(theta, y, w = rep(1, length(y))) {
  return(-sum(w * dnorm(log(y), mean = theta[1], sd = theta[2], log = TRUE)))
}



opt <- optim(c(0, 1), function(x) f(x, y, w), method = "L-BFGS-B", lower = c(-Inf, 1e-12), upper = c(Inf, Inf), hessian = TRUE)
diag(opt$hess)
opt <- optim(c(0, 1), function(x) f(x, y, w), 
             gr = function(x) rowSums(grads(x, y, w)),
             method = "L-BFGS-B", lower = c(-Inf, 1e-12), upper = c(Inf, Inf), 
             hessian = TRUE)

sqrt(diag(solve(opt$hess)))
rowSums(grads(opt$par, y, w))
rowSums(grads(theta, y, w))

diag(hess(theta, y, w))

vcov(lm(log(y)~1, weights = w))



# WR-design is almost E-optimal for WOR-design.
res_wor <- calculate_sampling_scheme(rep(1, N), n, G, H, w, crit = "E", design = "PO-WOR", tol = 5)
res_wr <-  calculate_sampling_scheme(rep(1, N), n, G, H, w, crit = "E", design = "PO-WR")
plot(res_wr$mu, res_wor$mu, bty = "l")




# D- and E-optimality, start from random scheme.
res1 <- calculate_sampling_scheme(n, G, H, w, crit = "D", design = "PO-WOR", init = runif(N))
res2 <- calculate_sampling_scheme(n, G, H, w, crit = "D", design = "PO-WR" , init = runif(N))
res3 <- calculate_sampling_scheme(n, G, H, w, crit = "E", design = "PO-WOR", init = runif(N))
res4 <- calculate_sampling_scheme(n, G, H, w, crit = "E", design = "PO-WOR", init = runif(N), tol = 5)
res5 <- calculate_sampling_scheme(n, G, H, w, crit = "E", design = "PO-WR" , init = runif(N))



res7 <- calculate_sampling_scheme(n, G_orth2,     H_orth2,     w, crit = "D")
res8 <- calculate_sampling_scheme(n, G_trans2,    H_trans2,    w, crit = "D")



# Orthogonalised with respect to gradient SSCP matrix. 
GtG <- matrix(0, nrow(G), nrow(G))
for ( i in 1:ncol(G) ) {
  GtG <- GtG + tcrossprod(G[, i])
}
A <- sqrtm(GtG)
Y_orth2 <- t(solve(A, t(Y)))
theta_orth2 <- colSums(w * Y_orth2) / sum(w)
G_orth2 <- grads(theta_orth2, Y_orth2, w)
H_orth2 <- hess(theta_orth2, w)
theta_orth2
t(solve(A) %*% theta)


# Orthogonalised with respect to gradient SSCP matrix + coordinate-wise affine transformation.
a <- rnorm(3)
b <- rnorm(3)
Y_trans2 <- t(a + t(Y_orth2) * b)
theta_trans2 <- colSums(w * Y_trans2) / sum(w)
G_trans2 <- grads(theta_trans2, Y_trans2, w)
H_trans2 <- hess(theta_trans2, w)
theta_trans2
a + b * theta_orth2



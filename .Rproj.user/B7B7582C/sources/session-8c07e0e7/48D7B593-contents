################################################################################
#
# File name: main.R
#
# Author: Henrik Imberg
#
# Last edited: 2023-03-17
#
# Description: Main function for the experiments in Section 6.
# 
# INPUT: Virtual simulation data stored in the Data folder.
#
# OUTPUT: is stored in the Tables folder.
#
################################################################################

# Clean up. ----
rm(list = ls())
cat("\14")


# Load packages. ----

library("expm")
library("sampling")
library("tidyverse")
library("xtable")

source("Rscript/acov.R")
source("Rscript/fmt.R")
source("Rscript/Lopt.R")
source("Rscript/Phi.R")
source("Rscript/Phiopt.R")


# Prepare data. ----

load("Data/VirtSim.R")
df$caseID <- as.numeric(df$caseID)

maximpact <- df %>% 
  group_by(caseID) %>% 
  summarise(max_impact_speed = max(impact_speed0, na.rm = TRUE), .groups = "keep") %>% 
  ungroup() %>% 
  dplyr::select(caseID, max_impact_speed)

VSdata <- df %>%
  dplyr::rename("OEOFF" = eoff,
                "prob" = eoff_acc_prob) %>%
  mutate(dec = -acc, 
         crash0 = as.numeric(impact_speed0 > 0),
         crash1 = as.numeric(impact_speed1 > 0),
         impact_speed_reduction = impact_speed0 - impact_speed1,
         injury_risk_reduction = injury_risk0 - injury_risk1,
         crash_avoidance = (1 - crash1) * crash0) %>%
  dplyr::select(caseID, OEOFF, dec, everything(), -acc) %>% 
  left_join(maximpact, by = "caseID")



# Baseline impact speed distribution. ----
cat("\n------------------ Baseline impact speed distribution ------------------\n\n")

# Clean up.
rm(list = setdiff(ls(), c("VSdata", lsf.str())))

# Data.
data <- VSdata %>% 
  filter(crash0 == 1)

# As vector. 
y <- data$impact_speed0 
logy <- log(data$impact_speed0) 

N <- nrow(data) # Full-data size.
n <- round(N / 100) # Subsample size.
w <- with(data, crash0 * prob) # Observation weights.

# Full-data parameter.
theta0 <- c(weighted.mean(logy, w), sqrt(cov.wt(matrix(logy, ncol = 1), w, method = "ML")$cov)) 

# Gradients. 
grads <- function(theta0, y, w = rep(1, nrow(y))) {
  Z <- (log(y) - theta0[1]) / theta0[2]
  dldMu <- - w * Z / theta0[2]
  dldSig <- - w * (Z^2 - 1) / theta0[2]
  return(rbind(dldMu, dldSig))
}
G <- grads(theta0, y, w)

# Hessian.
hess <- function(theta0, y, w) {
  d1 <- sum(w)
  return(diag(sum(w) * c(1, 2)) / theta0[2]^2)
}
H <- hess(theta0, y, w)

# V-matrix.
V <- matrix(0, nrow(H), ncol(H))
for ( i in 1:nrow(G) ) {
  V <- V + tcrossprod(G[, i])
}

# Parameters.
params <- tibble(crit = c("A", rep("c", 2), "D", rep("L", 2), "E", rep("Phi_r", 3)), 
                 name = c("A", "c, $\\mathbf{c} = (1,0)\\T$", "c, $\\mathbf{c} = (0,1)\\T$", "D", "d$^*_{\\mathrm{ER}}$", "d$^*_{\\mathrm{S}}$", "E", "$\\Phi_{0.5}$", "$\\Phi_5$", "$\\Phi_{10}$"), 
                 r = c(rep(NA, 7), 0.5, 5, 10),
                 design = "PO-WR")
L <- list(diag(ncol(H)), c(1, 0), c(0, 1), NULL, sqrtm(H), H %*% solve(sqrtm(V)), NULL, NULL, NULL, NULL)
params$L <- L
params$s <- 1:nrow(params)

# To store results. 
res <- add_column(params, niter = NA_real_, t = NA_real_, status = NA_real_,
                  A_eff = NA_real_, c1_eff = NA_real_, c2_eff = NA_real_, 
                  D_eff = NA_real_, dER_eff = NA_real_, dS_eff = NA_real_, E_eff = NA_real_, 
                  Phi05_eff = NA_real_, Phi5_eff = NA_real_, Phi10_eff = NA_real_)

# Find optimal sampling schemes and evaluate performance analytically. 
nreps <- 100
for ( i in 1:nrow(params) ) {
  t <- 0
  for ( j in 1:nreps ) { # Repeat to assess computation time.
    print_level <- ifelse( j == 1, 1, 0)
    if ( params$crit[i] %in% c("A", "c", "L") ) {
      opt <- Lopt(n, G, H, L = params$L[[i]], design = params$design[i], print_level = print_level)
      opt$status <- 0
      Gamma <- acov(mu = opt$mu, G, H, params$design[i])
    } else {
      opt <- Phiopt(n, G, H, crit = params$crit[i], 
                    L = params$L[[i]], r = params$r[i], design = params$design[i], 
                    print_level = print_level)
      Gamma <- opt$Gamma
    }
    t <- t + opt$elapsed / nreps
  }
  res$niter[i] <- opt$niter
  res$t[i] <- t
  res$status[i] <- opt$status
  res$A_eff[i] <- Phi(Gamma, crit = "A")
  res$c1_eff[i] <- Phi(Gamma, crit = "c", L = L[[2]])
  res$c2_eff[i] <- Phi(Gamma, crit = "c", L = L[[3]])
  res$D_eff[i] <- Phi(Gamma, crit = "D")
  res$dER_eff[i] <- Phi(Gamma, crit = "L", L = L[[5]])
  res$dS_eff[i] <- Phi(Gamma, crit = "L", L = L[[6]])
  res$E_eff[i] <- Phi(Gamma, crit = "E")
  res$Phi05_eff[i] <- Phi(Gamma, crit = "Phi_r", r = 0.5)
  res$Phi5_eff[i] <- Phi(Gamma, crit = "Phi_r", r = 5)
  res$Phi10_eff[i] <- Phi(Gamma, crit = "Phi_r", r = 10)
}

# Relative efficiencies.
res <- res %>% 
  mutate(A_eff = min(A_eff) / A_eff, 
         c1_eff = min(c1_eff) / c1_eff, 
         c2_eff = min(c2_eff) / c2_eff, 
         D_eff = min(D_eff) / D_eff,
         E_eff = min(E_eff) / E_eff,
         dER_eff = min(dER_eff) / dER_eff,
         dS_eff = min(dS_eff) / dS_eff, 
         Phi05_eff = min(Phi05_eff) / Phi05_eff,
         Phi5_eff = min(Phi5_eff) / Phi5_eff,
         Phi10_eff = min(Phi10_eff) / Phi10_eff) %>% 
  mutate(niter_s = ifelse(status == 1, "Diverged", sprintf("%.0f", niter)),
         niter_s = ifelse(status == 2, "Did not converge", niter_s),
         niter_s = ifelse(niter_s == "NA", "0", niter_s)) %>% 
  dplyr::select(crit, name, niter, niter_s, t, everything())

# If optimal sampling scheme could not be found.
for ( i in 1:nrow(res) ) {
  if ( res$status[i] > 0 ) {
    res[i, 11:ncol(res)] <- NA
    if ( res$crit[i] == "E" ) {
      res$E_eff[which(res$design == res$design[i])] <- NA
    } else if ( res$crit[i] == "Phi_r" & res$r[i] == 5 ) {
      res$Phi5_eff[which(res$design == res$design[i])] <- NA
    } else if ( res$crit[i] == "Phi_r" & res$r[i] == 10 ) {
      res$Phi10_eff[which(res$design == res$design[i])] <- NA
    }
  }
}

# Format results and write as xtable to file.
res %>% mutate(across(contains("eff"), fmt)) %>% 
  mutate(t = ifelse(status == 0, sprintf("%.2f", t), "-")) %>% 
  dplyr::select(name, niter_s, t, A_eff, c1_eff, c2_eff, D_eff, dER_eff, Phi5_eff) %>% 
  dplyr::rename("Optimality criterion." = name,
                "No. iterations" = niter_s,
                "Time (s)" = t,
                "A-eff" = A_eff,
                "c$_{(1,0)}$-eff" = c1_eff,
                "c$_{(0,1)}$-eff" = c2_eff,
                "D-eff" = D_eff,
                "d$^*_{\\mathrm{ER}}$-eff"  = dER_eff,
                "$\\Phi_5$-eff" = Phi5_eff) %>% 
  xtable(caption = "Number of iterations needed for convergence, average execution time, and relative efficiencies of sampling schemes and optimality criteria for estimating the log-normal model \\eqref{eq:lognormal}. \\textit{eff = relative efficiency.}", 
         label = "tab:baseline_impact_speed", 
         align = c("llrrrrrrrr")) %>% 
  print(digits = 2, 
        caption.placement = "top", 
        hline.after = c(0, nrow(.)),
        include.rownames = FALSE,
        sanitize.text.function = function(x) x,
        table.placement = "htb!", 
        file = "Tables/tab2.tex")

# Execution time relatively to D-optimality.
cat("Relative execution time of L-optimality compared to D-optimality.\n")
res %>% 
  mutate(t_rel = sprintf("%.2f", t / t[which(crit == "D")])) %>% 
  filter(crit %in% c("A", "c", "L")) %>% 
  dplyr::select(name, t_rel) %>% 
  print(decimals = 0)



# Baseline impact speed response surface. ----
cat("\n------------------ Baseline impact response surface ------------------\n\n")

# Clean up.
rm(list = setdiff(ls(), c("VSdata", lsf.str())))

# Data.
data <- VSdata 
data$caseID <- factor(data$caseID)

N <- nrow(data) # Full-data size.
n <- round(N / 100) # Subsample size.

# Model matrix and response vector.
X <- model.matrix(~caseID*dec*OEOFF, data = data)
y <- with(data, impact_speed0 / max_impact_speed)

# Fit quasi-binomial logistic response model.
t0 <- Sys.time()
fit <- glm(y~caseID*dec*OEOFF, data = data, family = "quasibinomial")
t1 <- Sys.time()
t_theta0 <- t1 - t0

theta0 <- coef(fit) # Full-data parameter.
p <- as.numeric(predict(fit, type = "response")) # Predictions.
G <- t((y - p) * X) # Gradients
H <- solve(vcov(fit)) # Hessian

# V-matrix.
V <- matrix(0, nrow(H), ncol(H))
for ( i in 1:nrow(X) ) {
  V <- V + tcrossprod(G[, i])
}

# Parameters.
params <- tibble(crit = c("A", "D", rep("L", 2), "E", rep("Phi_r", 3)), 
                 name = c("A", "D", "d$^*_{\\mathrm{ER}}$", "d$^*_{\\mathrm{S}}$", "E", "$\\Phi_{0.5}$", "$\\Phi_5$", "$\\Phi_{10}$"), 
                 r = c(rep(NA, 5), 0.5, 5, 10),
                 design = "PO-WR")
L <- list(diag(ncol(H)), NULL, sqrtm(H), H %*% solve(sqrtm(V)), NULL, NULL, NULL, NULL)
params$L <- L
params$s <- 1:nrow(params)

# To store results.
res <- add_column(params, niter = NA_real_, t = NA_real_, status = NA_real_, 
                  A_eff = NA_real_, D_eff = NA_real_, 
                  dER_eff = NA_real_, dS_eff = NA_real_, E_eff = NA_real_, 
                  Phi05_eff = NA_real_, Phi5_eff = NA_real_, Phi10_eff = NA_real_)

# Find optimal sampling schemes and evaluate performance analytically. 
nreps <- 100
for ( i in 1:nrow(params) ) {
  t <- 0
  for ( j in 1:nreps ) { # Repeat to assess computation time.
    print_level <- ifelse( j == 1, 1, 0)
    if ( params$crit[i] %in% c("A", "c", "L") ) {
      opt <- Lopt(n, G, H, L = params$L[[i]], design = params$design[i], print_level = print_level)
      opt$status <- 0
      Gamma <- acov(mu = opt$mu, G, H, params$design[i])
    } else {
      opt <- Phiopt(n, G, H, crit = params$crit[i], L = params$L[[i]], r = params$r[i], design = params$design[i], print_level = print_level)
      Gamma <- opt$Gamma
    }
    t <- t + opt$elapsed / nreps
  }
  res$niter[i] <- opt$niter
  res$t[i] <- t
  res$status[i] <- opt$status
  res$A_eff[i] <- Phi(Gamma, crit = "A")
  res$D_eff[i] <- Phi(Gamma, crit = "D")
  res$dER_eff[i] <- Phi(Gamma, crit = "L", L = L[[3]])
  res$dS_eff[i] <- Phi(Gamma, crit = "L", L = L[[4]])
  res$E_eff[i] <- Phi(Gamma, crit = "E")
  res$Phi05_eff[i] <- Phi(Gamma, crit = "Phi_r", r = 0.5)
  res$Phi5_eff[i] <- Phi(Gamma, crit = "Phi_r", r = 5)
  res$Phi10_eff[i] <- Phi(Gamma, crit = "Phi_r", r = 10)
}

# Relative efficiencies.
res <- res %>% 
  mutate(A_eff = min(A_eff) / A_eff, 
         D_eff = min(D_eff) / D_eff,
         E_eff = min(E_eff) / E_eff,
         dER_eff = min(dER_eff) / dER_eff,
         dS_eff = min(dS_eff) / dS_eff, 
         Phi05_eff = min(Phi05_eff) / Phi05_eff,
         Phi5_eff = min(Phi5_eff) / Phi5_eff,
         Phi10_eff = min(Phi10_eff) / Phi10_eff) %>% 
  mutate(niter_s = ifelse(status == 1, "Diverged", sprintf("%.0f", niter)),
         niter_s = ifelse(status == 2, "Did not converge", niter_s),
         niter_s = ifelse(niter_s == "NA", "0", niter_s)) %>% 
  dplyr::select(crit, name, niter, niter_s, t, everything())

# If optimal sampling scheme could not be found.
for ( i in 1:nrow(res) ) {
  if ( res$niter_s[i] == "Did not converge") {
    res[i, 11:ncol(res)] <- NA
    if ( res$crit[i] == "E" ) {
      res$E_eff[which(res$design == res$design[i])] <- NA
    } else if ( res$crit[i] == "Phi_r" & res$r[i] == 5 ) {
      res$Phi5_eff[which(res$design == res$design[i])] <- NA
    } else if ( res$crit[i] == "Phi_r" & res$r[i] == 10 ) {
      res$Phi10_eff[which(res$design == res$design[i])] <- NA
    }
  }
}

# Format results and write as xtable to file.
res %>% mutate(across(contains("eff"), fmt)) %>% 
  mutate(t = ifelse(status == 0, sprintf("%.2f", t), "-")) %>% 
  dplyr::select(name, niter_s, t, A_eff, D_eff, dER_eff, dS_eff, Phi05_eff) %>% 
  dplyr::rename("Optimality criterion." = name,
                "No. iterations" = niter_s,
                "Time (s)" = t,
                "A-eff" = A_eff,
                "D-eff" = D_eff,
                "d$^*_{\\mathrm{ER}}$-eff"  = dER_eff,
                "d$^*_{\\mathrm{S}}$-eff"  = dS_eff,
                "$\\Phi_{0.5}$-eff" = Phi05_eff) %>% 
  xtable(caption = sprintf("Number of iterations needed for convergence, average execution time, and relative efficiencies of sampling schemes and optimality criteria for estimating the quasi-binomial logistic regression model \\eqref{eq:qblr}. The computation time for fitting the model to the full dataset was %.2f s. \\textit{eff = relative efficiency.}", t_theta0), 
         label = "tab:impact_speed_response_surface", 
         align = c("llrrrrrrr")) %>% 
  print(digits = 2, 
        caption.placement = "top", 
        hline.after = c(0, nrow(.)),
        include.rownames = FALSE,
        sanitize.text.function = function(x) x,
        table.placement = "htb!", 
        file = "Tables/tab3.tex")

# Execution time relatively to D-optimality.
cat("Relative execution time of L-optimality compared to D-optimality.\n")
res %>% 
  mutate(t_rel = sprintf("%.2f", t / t[which(crit == "D")])) %>% 
  filter(crit %in% c("A", "c", "L")) %>% 
  dplyr::select(name, t_rel) %>% 
  print(decimals = 0)


# AEB vs. manual baseline driving. ----
cat("\n------------------ AEB vs. manual baseline driving ------------------n\n")

# Clean up.
rm(list = setdiff(ls(), c("VSdata", lsf.str())))

# Data.
data <- VSdata %>% 
  filter(crash0 == 1) 

# As matrix.
Y <- data %>% 
  dplyr::select(impact_speed_reduction, injury_risk_reduction, crash_avoidance) %>%
  as.matrix()

N <- nrow(data) # Full-data size.
n <- round(N / 100) # Subsample size.
w <- with(data, crash0 * prob) # Observation weights.

theta0 <- colSums(w * Y) / sum(w) # Full-data parameter.
G <- t(-w * t(t(Y) - theta0)) # Gradients.
H <- diag(sum(w), length(theta0)) # Hessian.

# V-matrix.
V <- matrix(0, nrow(H), ncol(H))
for ( i in 1:ncol(G) ) {
  V <- V + tcrossprod(G[, i])
}

# Parameters.
params <- tibble(crit = c("A", rep("c", 3), "D", rep("L", 2), "E", rep("Phi_r", 3)), 
                 name = c("A", "c, $\\mathbf{c} = (1,0,0)\\T$", "c, $\\mathbf{c} = (0,1,0)\\T$", "c, $\\mathbf{c} = (0,0,1)\\T$", "D", "d$^*_{\\mathrm{ER}}$", "d$^*_{\\mathrm{S}}$", "E", "$\\Phi_{0.5}$", "$\\Phi_5$", "$\\Phi_{10}$"), 
                 r = c(rep(NA, 8), 0.5, 5, 10),
                 design = "PO-WR")
L <- list(diag(ncol(H)), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), NULL, sqrtm(H), H %*% solve(sqrtm(V)), NULL, NULL, NULL, NULL)
params$L <- L
params$s <- 1:nrow(params)

# To store results.
res <- add_column(params, niter = NA_real_, t = NA_real_, status = NA_real_, 
                  A_eff = NA_real_, c1_eff = NA_real_, c2_eff = NA_real_, c3_eff = NA_real_, 
                  D_eff = NA_real_, dER_eff = NA_real_, dS_eff = NA_real_, E_eff = NA_real_, 
                  Phi05_eff = NA_real_, Phi5_eff = NA_real_, Phi10_eff = NA_real_)

# Find optimal sampling schemes and evaluate performance analytically. 
nreps <- 10000
for ( i in 1:nrow(params) ) {
  t <- 0
  for ( j in 1:nreps ) { # Repeat to assess computation time.
    print_level <- ifelse( j == 1, 1, 0)
    if ( params$crit[i] %in% c("A", "c", "L") ) {
      opt <- Lopt(n, G, H, L = params$L[[i]], design = params$design[i], print_level = print_level)
      opt$status <- 0
      Gamma <- acov(mu = opt$mu, G, H, params$design[i])
    } else {
      opt <- Phiopt(n, G, H, crit = params$crit[i], L = params$L[[i]], r = params$r[i], design = params$design[i], print_level = print_level)
      Gamma <- opt$Gamma
    }
    t <- t + opt$elapsed / nreps
  }
  res$niter[i] <- opt$niter
  res$t[i] <- t
  res$status[i] <- opt$status
  res$A_eff[i] <- Phi(Gamma, crit = "A")
  res$c1_eff[i] <- Phi(Gamma, crit = "c", L = L[[2]])
  res$c2_eff[i] <- Phi(Gamma, crit = "c", L = L[[3]])
  res$c3_eff[i] <- Phi(Gamma, crit = "c", L = L[[4]])
  res$D_eff[i] <- Phi(Gamma, crit = "D")
  res$dER_eff[i] <- Phi(Gamma, crit = "L", L = L[[6]])
  res$dS_eff[i] <- Phi(Gamma, crit = "L", L = L[[7]])
  res$E_eff[i] <- Phi(Gamma, crit = "E")
  res$Phi05_eff[i] <- Phi(Gamma, crit = "Phi_r", r = 0.5)
  res$Phi5_eff[i] <- Phi(Gamma, crit = "Phi_r", r = 5)
  res$Phi10_eff[i] <- Phi(Gamma, crit = "Phi_r", r = 10)
}

# Relative efficiencies.
res <- res %>% 
  mutate(A_eff = min(A_eff) / A_eff, 
         c1_eff = min(c1_eff) / c1_eff, 
         c2_eff = min(c2_eff) / c2_eff, 
         c3_eff = min(c3_eff) / c3_eff, 
         D_eff = min(D_eff) / D_eff,
         E_eff = min(E_eff) / E_eff,
         dER_eff = min(dER_eff) / dER_eff,
         dS_eff = min(dS_eff) / dS_eff, 
         Phi05_eff = min(Phi05_eff) / Phi05_eff,
         Phi5_eff = min(Phi5_eff) / Phi5_eff,
         Phi10_eff = min(Phi10_eff) / Phi10_eff) %>% 
  mutate(niter_s = ifelse(status == 1, "Diverged", sprintf("%.0f", niter)),
         niter_s = ifelse(status == 2, "Did not converge", niter_s),
         niter_s = ifelse(niter_s == "NA", "0", niter_s)) %>% 
  dplyr::select(crit, name, niter, niter_s, t, everything())

# If optimal sampling scheme could not be found.
for ( i in 1:nrow(res) ) {
  if ( res$niter_s[i] == "Did not converge") {
    res[i, 11:ncol(res)] <- NA
    if ( res$crit[i] == "E" ) {
      res$E_eff[which(res$design == res$design[i])] <- NA
    } else if ( res$crit[i] == "Phi_r" & res$r[i] == 5 ) {
      res$Phi5_eff[which(res$design == res$design[i])] <- NA
    } else if ( res$crit[i] == "Phi_r" & res$r[i] == 10 ) {
      res$Phi10_eff[which(res$design == res$design[i])] <- NA
    }
  }
}

# Format results and write as xtable to file.
res %>% mutate(across(contains("eff"), fmt)) %>% 
  mutate(t = ifelse(status == 0, sprintf("%.2f", t), "-")) %>% 
  dplyr::select(name, niter_s, t, A_eff, c1_eff, c2_eff, c3_eff, D_eff, E_eff) %>% 
  dplyr::rename("Optimality criterion." = name,
                "No. iterations" = niter_s,
                "Time (s)" = t,
                "A-eff" = A_eff,
                "c$_{(1,0,0)}$-eff" = c1_eff,
                "c$_{(0,1,0)}$-eff" = c2_eff,
                "c$_{(0,0,1)}$-eff" = c3_eff,
                "D-eff" = D_eff,
                "E-eff" = E_eff) %>% 
  xtable(caption = "Number of iterations needed for convergence, average execution time, and relative efficiencies of sampling schemes and optimality criteria for estimating the vector of finite population means \\eqref{eq:theta0_fps}. \\textit{eff = relative efficiency.}", 
         label = "tab:finite_population_inference", 
         align = c("llrrrrrrrr")) %>% 
  print(digits = 2,
        caption.placement = "top", 
        hline.after = c(0, nrow(.)),
        include.rownames = FALSE,
        sanitize.text.function = function(x) x,
        table.placement = "htb!", 
        file = "Tables/tab4.tex")

# Execution time relatively to D-optimality.
cat("Relative execution time of L-optimality compared to D-optimality.\n")
res %>% 
  mutate(t_rel = sprintf("%.2f", t / t[which(crit == "D")])) %>% 
  filter(crit %in% c("A", "c", "L")) %>% 
  dplyr::select(name, t_rel) %>% 
  print(decimals = 0)

rm(list = ls())
cat("\14")

library("nloptr")
library("scales") # to access break formatting functions
library("tictoc")
library("tidyverse")

set.seed(123)

Nseq <- 10^seq(1, 5, length.out = 21)
yseq <- rexp(max(Nseq))

niter <- status <- f_rel_diff <- rep(NA, length = length(Nseq))
tictoc <- matrix(NA, nrow =length(Nseq), ncol = 2)
colnames(tictoc) <- c("Analytical", "Numerical")

# NLopt options.
opts <- list("algorithm" = "NLOPT_LD_AUGLAG_EQ",
             "maxeval" = 1e3, 
             "print_level" = 0,
             "local_opts" = list("algorithm" = "NLOPT_LD_MMA", 
                                 "maxeval" = 1e3))

for ( i in seq_along(Nseq) ) {
  
  N <- Nseq[i]
  y <- yseq[1:N]
  n <- floor(max(1, N / 10))
  
  # Analytical.
  tic()
  mu <- n * y / sum(y)
  t <- toc(quiet = TRUE)
  tictoc[i, 1] <- as.numeric(t$toc - t$tic)
  
  # Numerical.
  tic()
  opt <- nloptr( x0 = rep(n / N, N),
                 eval_f = function(mu) sum(y^2 / mu),
                 eval_grad_f = function(mu) -y^2 / mu^2,
                 lb = rep(0, N), 
                 ub = rep(n, N), 
                 eval_g_eq = function(mu) sum(mu) - n, 
                 eval_jac_g_eq = function(mu) rep(1, length(mu)), 
                 opts = opts)
  t <- toc(quiet = TRUE)
  status[i] <- opt$status
  niter[i] <- opt$iterations
  tictoc[i, 2] <- as.numeric(t$toc - t$tic)
  
}

res <- tibble(N = Nseq,
              niter = niter,
              status = status,
              tictoc_analytical = tictoc[, 1], 
              tictoc_numerical = tictoc[, 2])

# ggplot(data = res) + 
#   geom_line(aes(x = N, y = tictoc_analytical, colour = "Analytical", lty = "Analytical")) +   
#   geom_line(aes(x = N, y = tictoc_numerical, colour = "Numerical", lty = "Numerical")) +   
#   geom_point(aes(x = N, y = tictoc_analytical, colour = "Analytical")) +   
#   geom_point(data = res %>% filter(status == 4), aes(x = N, y = tictoc_numerical, colour = "Numerical")) + 
#   geom_point(data = res %>% filter(status == 5), aes(x = N, y = tictoc_numerical, colour = "Numerical", shape = as.character(status))) + 
#   scale_color_brewer(palette = "Dark2") + 
#   scale_x_continuous(trans = "log10") +
#   scale_y_continuous(breaks = c(seq(0, 240, 60), seq(300, 1500, 300)), labels = c("0 sec", rep("", 4), "5 min", "10 min", "15 min", "20 min", "25 min")) + 
#   scale_linetype_manual(values = c(2, 1)) +
#   scale_shape_manual(values = 4, labels = "Did not converge") +
#   labs(x = "Dataset size",
#        y = "Time", 
#        colour = NULL,
#        lty = NULL,
#        shape = NULL) + 
#   guides(colour = guide_legend(order = 1), 
#          lty = guide_legend(order = 1),
#          shape = guide_legend(order = 2)) + 
#   theme(legend.justification = c(0, 1), 
#         legend.position = c(0, 1))

ptsize <- 10
theme_set(theme_classic()) 
theme_update(axis.text = element_text(size = ptsize, colour = "black", family = "serif"),
             axis.line = element_line(colour = "black", linewidth = 0.25), 
             axis.ticks = element_line(colour = "black", linewidth = 0.25), 
             axis.title.x = element_text(margin = margin(t = 0.25, r = 0, b = 0, l = 0, unit = 'cm')),
             axis.title.y = element_text(margin = margin(t = 0, r = 0.25, b = 0, l = 0, unit = 'cm')),
             legend.background = element_blank(),
             legend.key.width = unit(0, "cm"),
             legend.key.height = unit(0.4, "cm"),
             legend.margin = margin(0, 0, 0, 0.5, "cm"),
             legend.spacing =  unit(0, "cm"),
             legend.position = "bottom",
             legend.text = element_text(size = ptsize, colour = "black", family = "serif"),
             legend.title = element_text(size = ptsize, colour = "black", family = "serif"),
             strip.background.x = element_blank(),
             panel.border = element_blank(),
             panel.grid = element_blank(),  
             plot.subtitle = element_text(size = ptsize, colour = "black", family = "serif", face = "plain", hjust = 0),
             plot.title = element_text(size = ptsize, colour = "black", family = "serif", face = "plain", hjust = 0),
             text = element_text(size = ptsize, colour = "black", family = "serif"))

ggplot(data = res, aes(x = N, y = tictoc_numerical)) + 
  geom_line() +   
  geom_point(data = res %>% filter(status == 4), colour = "white", size = 2) + 
  geom_point(data = res %>% filter(status == 4)) + 
  geom_point(data = res %>% filter(status == 5), colour = "white", size = 3.5) + 
  geom_point(data = res %>% filter(status == 5), aes(shape = as.character(status))) + 
  scale_x_continuous(trans = "log10", labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(breaks = c(seq(0, 20, 5))) + 
  scale_shape_manual(values = 4, labels = "Did not converge") +
  labs(x = "Dataset size (N)",
       y = "Computation time (s)", 
       colour = NULL,
       lty = NULL,
       shape = NULL) + 
  guides(colour = guide_legend(order = 1), 
         lty = guide_legend(order = 1),
         shape = guide_legend(order = 2)) + 
  theme(legend.justification = c(0, 1), 
        legend.position = c(0, 1))

ggsave("Figures/nloptr.png", dpi = 1000, width = 80, height = 60, unit = "mm")

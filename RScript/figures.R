################################################################################
#
# File name: figures.R
#
# Author: Henrik Imberg
#
# Last edited: 2023-04-05
#
# Description: Generate Figure S1 and S2 in Appendix B.
# 
# INPUT: Virtual simulation data stored in the Data folder.
#
# OUTPUT: is stored in the Figures folder.
#
################################################################################

# Clean up. ----
rm(list = ls())
cat("\14")


# Load packages. ----

library("cowplot")
library("tidyverse")


# ggplot theme. ----

ptsize <- 10
theme_set(theme_classic()) 
theme_update(axis.text = element_text(size = ptsize, colour = "black", family = "serif"),
             axis.line = element_line(colour = "black", linewidth = 0.2), 
             axis.ticks = element_line(colour = "black", linewidth = 0.2), 
             axis.title.x = element_text(margin = margin(t = 0.25, r = 0, b = 0, l = 0, unit = 'cm')),
             axis.title.y = element_text(margin = margin(t = 0, r = 0.25, b = 0, l = 0, unit = 'cm')),
             legend.background = element_blank(),
             legend.key.width = unit(1, "cm"),
             legend.key.height = unit(0.4, "cm"),
             legend.margin = margin(0, 0, 0, 0.5, "cm"),
             legend.spacing =  unit(0, "cm"),
             legend.position = "bottom",
             legend.text = element_text(size = ptsize, colour = "black", family = "serif"),
             legend.title = element_text(size = ptsize, colour = "black", family = "serif"),
             strip.background = element_blank(),
             strip.text = element_blank(),
             panel.border = element_blank(),
             panel.grid = element_blank(),  
             plot.subtitle = element_text(size = ptsize, colour = "black", family = "serif", face = "plain", hjust = 0),
             plot.title = element_text(size = ptsize, colour = "black", family = "serif", face = "plain", hjust = 0),
             text = element_text(size = ptsize, colour = "black", family = "serif"))
update_geom_defaults("text", list(size = ptsize / ggplot2:::.pt))


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

# Data.
data <- VSdata %>% 
  filter(crash0 == 1)
w <- with(data, crash0 * prob) # Observation weights.


# Figure S1. 
figS1a <- ggplot(data, aes(x = impact_speed0)) + 
  geom_histogram(aes(y = after_stat(density), weight = w), bins = 30, colour = "black", fill = "grey70", linewidth = 0.2) + 
  labs(x = "Baseline impact speed (km/h)", 
       y = "Density")

figS1b <- ggplot(data, aes(x = impact_speed_reduction, y = after_stat(density), weight = w)) + 
  geom_histogram(bins = 30, colour = "black", fill = "grey70", linewidth = 0.2) + 
  labs(x = "Impact speed reduction (km/h)", 
       y = "Density")

figS1c <- ggplot(data, aes(x = injury_risk_reduction, weight = w)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 100, colour = "black", fill = "grey70", linewidth = 0.2) + 
  geom_rug(length = unit(0.02, "npc"), linewidth = 0.2) + 
  labs(x = "Injury risk reduction", 
       y = "Density")

figS1d <- ggplot(data, aes(x = crash_avoidance, weight = w)) + 
  geom_bar(aes(y = after_stat(prop)), width = 0.5, colour = "black", fill = "grey70", linewidth = 0.2) + 
  scale_x_continuous(breaks = c(0, 1), labels = c("0 (No)", "1 (Yes)")) + 
  scale_y_continuous(labels = scales::percent) + 
  labs(x = "Crash avoidance",
       y = "Proportion")

plot_grid(figS1a, figS1b, figS1c, figS1d, nrow = 2, ncol = 2, labels = LETTERS[1:4], label_fontfamily = "serif", label_size = ptsize)

ggsave("Figures/FigureS1.png", dpi = 1000, width = 159, height = 85, unit = "mm")



# Baseline impact speed response surface. ----

# Data.
data <- VSdata %>% 
  filter(caseID <= 10)
data$caseID <- factor(data$caseID)

# Fit quasi-binomial logistic response model.
y <- with(data, impact_speed0 / max_impact_speed)
fit <- glm(y~caseID*dec*OEOFF, data = data, family = "quasibinomial")
p <- as.numeric(predict(fit, type = "response")) # Predictions.

# Figure S2.
plt <- data %>% 
  mutate(impact_speed  = max_impact_speed * as.numeric(p), grp = 1) %>%
  bind_rows(data %>% mutate(impact_speed = impact_speed0, grp = 0)) %>% 
  filter(caseID %in% c(1, 2, 10)) 

figS2 <- ggplot(plt) +
  geom_line(aes(x = OEOFF, y = impact_speed, group = dec, colour = dec))  +
  geom_text(data = data.frame(grp = 0, caseID = as.factor(1), x = 0, y = 70, label = "Observed"), aes(x = x, y = y, label = label), hjust = 0, family = "serif") + 
  geom_text(data = data.frame(grp = 1, caseID = as.factor(1), x = 0, y = 70, label = "Predicted"), aes(x = x, y = y, label = label), hjust = 0, family = "serif") + 
  labs(x = "Off-road glance duration (s)",
       y = "Impact speed (km/h)",
       colour = expression(paste("Maximal deceleration (m/", s^2, ")"))) + 
  facet_grid(grp~caseID, labeller = labeller(grp = c("0" = "Observed", "1" = "Predicted"))) +
  theme(legend.position = "bottom") 

ggsave("Figures/FigureS2.png", dpi = 1000, width = 159, height = 85, unit = "mm")
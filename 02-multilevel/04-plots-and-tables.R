#####################################################
### Create plots and tables #########################
#####################################################


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### !!! Make sure that you are in the correct working directory! ####
###                     analysis/02-multilevel                   ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## The estimation needs to run first (there need to be logml.csv files in results/ml*/v_*/)
## If you pulled everything from Github, the logml files are already there. No need to estimate again to create the plots

library(tidyverse)

models = c("ml1" = "v_1", "ml2" = "v_1", "mlf" = "v_1") # If there are different versions for the different models, this can be changed here

nh = 8

helper = function(m) {
  tmp = readr::read_csv(sprintf("results/%s/%s/logml.csv", m, models[m]), col_names = FALSE)
  tmp = tmp$X1
  if (length(tmp) == (nh -1)) tmp = c(0, tmp)
  tmp
}

ml = sapply(names(models), helper)
ml[1,2:3] = ml[1,1] # Add the estimate for con1 to all models (its estimated in ml1)

plotdata = data.frame(ml = c(ml), hyp = rep(1:nh, times = 3), model = rep(names(models), each = nh))

lines_step = 50

# Plot 1

P = ggplot(plotdata, aes(x = hyp, y = ml, color = model, group = model)) +
  geom_line(linewidth = 0.5) +
  scale_x_continuous(name = "Hypothesis", breaks = 1:nh) +
  scale_y_continuous(name = "Log Marginal Likelihood",
                     limits = c(floor(min(ml)/lines_step)*lines_step, ceiling(max(ml)/lines_step)*lines_step),
                     breaks = seq(floor(min(ml)/lines_step)*lines_step, ceiling(max(ml)/lines_step)*lines_step, lines_step)) +
  scale_color_manual(name = "Model",
                     labels = c("mean con.", "full con.", "free"),
                     values = c("ml1" = "#9295ff",
                                "ml2" = "firebrick1",
                                "mlf" = "#000000")) +
  guides(color = guide_legend(position = "inside")) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "black", linewidth = 0.1),
        legend.position.inside = c(0.70,0.25),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))
P
#ggsave("ml-plot.png", P, dpi = 600, units = "cm", width = 8, height = 7)

# Plot 2

bf = plotdata$ml[plotdata$model== "ml1"] - plotdata$ml[plotdata$model== "mlf"]
plotdata2 = data.frame(bf = bf, hyp = 1:nh)
lines_step2 = 1

P2 = ggplot(plotdata2, aes(x = hyp, y = bf)) +
  geom_line(linewidth = 0.5) +
  scale_x_continuous(name = "Hypothesis", breaks = 1:nh) +
  scale_y_continuous(name = bquote(Log~BF~ML[mean~con]~vs.~ML[free]),
                     limits = c(floor(min(bf)/lines_step2)*lines_step2, ceiling(max(bf)/lines_step2)*lines_step2),
                     breaks = seq(floor(min(bf)/lines_step2)*lines_step2, ceiling(max(bf)/lines_step2)*lines_step2, lines_step2)) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "black", linewidth = 0.1),
        legend.position.inside = c(0.80,0.25),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))
P2
#ggsave("ml-plot2.png", P2, dpi = 600, units = "cm", width = 8, height = 7)

# Combined

Pc = gridExtra::grid.arrange(P, P2, nrow = 1)
#Pc = ggpubr::ggarrange(P, P2, nrow = 1, labels = c("A","B")) # for some reason, this makes the legend look bad. Added the labels with photoshop
ggsave("ml-plot-combined.png", Pc, dpi = 600, units = "cm", width = 16, height = 7)

## ------- Table

## ML
tdata = as.data.frame(t(ml))
tdata = tibble::add_column(tdata, Model = c("mean constraint", "full constraint", "free"), .before = 1)
colnames(tdata) = c("Model", sprintf("$\\mathcal{H}_%i$", 1:8))
rownames(tdata) = NULL

knitr::kable(tdata, digits = 2)

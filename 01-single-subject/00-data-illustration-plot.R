############################################
### Data Illustration plot (Figure 3) ######
############################################


library(tidyverse)

source("../get_data.R")

subject = 1

plotdata = filter(data, id == subject)

plotdata$shape = case_match(plotdata$reward,1 ~ 8, 0 ~ 1, -9 ~ NA)
plotdata$block = sprintf("Block %i", plotdata$pair)
P = ggplot(plotdata, aes(x = trial, y = choice, shape = as.factor(shape))) +
  geom_point() +
  facet_wrap(.~block, ncol = 2) +
  scale_shape_manual(values = c(1,8), name = "Outcome", labels = c("Negative", "Positive")) +
  scale_x_continuous(name = "Trial", minor_breaks = NULL) +
  scale_y_continuous(breaks = c(0,1), labels = c("B", "A"), limits = c(-1,2), name = "Choice",minor_breaks = NULL) +
  theme_classic()
P

ggsave("data-illustration-plot.png", P, width = 20, height = 10, units = "cm", dpi = 600)

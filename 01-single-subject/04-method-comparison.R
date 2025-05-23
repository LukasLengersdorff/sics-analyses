
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
##+++++++++++ SINGLE SUBJECT ANALYSES +++++++++++++++++++++++++++++++++++++++++
##+++++++++++     Method comparison   +++++++++++++++++++++++++++++++++++++++++
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### !!! Make sure that you are in the correct working directory! ####
###                     analysis/01-single-subject               ####
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


library(tidyverse)

# Load data
methods = c("uem", "cem", "sics")
cons = 1:8
datalist = list()

for (m in methods) {
  tmplist = list()
  for (con in cons) {
    tmplist[[sprintf("con%i", con)]] = safely(read_csv)(sprintf("results/%s/con%i.csv", m, con), show_col_types = FALSE)$result
  }
  datalist[[m]] = tmplist
}

# Retrieve p(D|C) from BFs for uem and cem

for (m in methods[1:2]) { #uem and cem are first two elements of `methods`
  for (con in names(datalist[[m]])) {
    datalist[[m]][[con]]$pDC = log(datalist[[m]][[con]]$bf) + datalist[[m]][[con]]$pD
  }
}

## Long data

longdata = list()

for (m in methods) {
  longdata[[m]] = list()
  for (con in names(datalist[[m]])) {
    longdata[[m]][[con]] = data.frame(method = rep(m, nrow(datalist[[m]][[con]])),
                                      con = rep(con, nrow(datalist[[m]][[con]])),
                                      ml = datalist[[m]][[con]]$pDC,
                                      t = datalist[[m]][[con]]$t.user + datalist[[m]][[con]]$t.system)
  }
  longdata[[m]] = list_rbind(longdata[[m]])
}
longdata = list_rbind(longdata)

longdata = mutate(longdata, method = factor(method, levels = c("uem", "cem", "sics")), hypothesis = str_replace(con, "con", ""))
longdata = mutate(longdata, hypothesis2 = sprintf("H[%s]",hypothesis))


# Just for the plot, we append copies of filter(longdata, method == "sics", con = "con1") with method changed to "cem" and "uem".
# This is because for hypothesis 1, all MLs were estimated solely by bridgesampling
for (m in c("uem", "cem")) {
  tmp = filter(longdata, method == "sics", con == "con1")
  tmp$method = m
  longdata = rbind(longdata, tmp)
}

plotdata = filter(longdata, is.finite(ml))

P = ggplot(data = plotdata, aes(x = method, y = ml, color = method)) +
  geom_violin(scale = "width") +
  geom_jitter(width = 0.05, alpha = 0.25, size = 0.2) +
  facet_wrap(.~hypothesis2, nrow = 1, scale = "fixed", labeller = label_parsed) +
  scale_x_discrete(name = "Estimation Method", labels = c("uEP", "cEP", "SICS")) +
  scale_y_continuous(name = "Log Marginal Likelihood",
                     #limits = c(floor(min(plotdata$ml)), ceiling(max(plotdata$ml))),
                     breaks = seq(floor(min(plotdata$ml)), ceiling(max(plotdata$ml)), 2)) +
  scale_color_manual(guide = NULL,
                       values = c("uem" = "#a5a8ff",
                                  "cem" = "#9295ff",
                                  "sics" = "#8184ff")) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "black", linewidth = 0.1),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background = element_rect(colour = NA, fill = "lightgrey"))
P

ggsave("comparison-plot.png", P, dpi = 600, width = 16.72, height = 6.14, units = "cm")


## Summary table

tdata = mutate(longdata, ml = if_else(is.finite(ml), ml, NaN))

agg = tdata %>% group_by(method, con) %>% summarize(method = unique(method), con = unique(con), m = mean(ml, na.rm = TRUE), sd = sd(ml, na.rm = TRUE), t = mean(t))
tab = pivot_wider(agg, id_cols = con, names_from = method, values_from = c(m, sd, t))
tab_print = tab[,c(1,2,5,8,3,6,9,4,7,10)]

tab_print$con = sprintf("$\\mathcal{H}_%i$", 1:8)
colnames(tab_print) = c("Hypothesis", "$M_{uEP}$", "$SD_{uEP}$", "$\\bar{t}_{uEP}$","$M_{cEP}$", "$SD_{cEP}$","$\\bar{t}_{cEP}$", "$M_{SICS}$", "$SD_{SICS}$", "$\\bar{t}_{SICS}$")

knitr::kable(tab_print, digits = c(0,1,3,1,1,3,1,1,3,1))

## Efficiency

tab_eff = with(tab, data.frame(con = con,
                               eff_uem = (1/sd_uem)/t_uem,
                               eff_cem = (1/sd_cem)/t_cem,
                               eff_sics = (1/sd_sics)/t_sics))
tab_eff = mutate(tab_eff, cem_v_uem = eff_cem/eff_uem,
                      sics_v_uem = eff_sics/eff_uem,
                      sics_v_cem = eff_sics/eff_cem)

tab_eff$con = sprintf("$\\mathcal{H}_%i$", 1:8)
colnames(tab_eff) = c("Hypothesis",
                      "$\\text{eff}_\\text{uEP}$",
                      "$\\text{eff}_\\text{cEP}$",
                      "$\\text{eff}_\\text{SICS}$",
                      "$\\text{eff}_\\text{cEP}/\\text{eff}_\\text{uEP}$",
                      "$\\text{eff}_\\text{SICS}/\\text{eff}_\\text{uEP}$",
                      "$\\text{eff}_\\text{SICS}/\\text{eff}_\\text{cEP}$")

knitr::kable(tab_eff, digits = c(0,3,3,3,1,1,1))

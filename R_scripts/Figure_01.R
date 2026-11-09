# ============================================================================
# Figure 1 -- R panels
#
# Run source("_setup.R") first for the shared library() calls.
# ============================================================================


# --------------------------------------------------------------------------------------------------------------
# Figure 1 panel a -- Clinical annotation heatmap (Diagnosis, IPSS-M, Age, IMPACT, SF3B1) across Xenium samples
# --------------------------------------------------------------------------------------------------------------

df <- read.delim('data/Fig1a.csv', sep = ',')

df$Age.cat <- '40-60'
df$Age.cat[which(df$Age < 40)] <-'Under 40'
df$Age.cat[which(df$Age > 60)] <- 'Over 60'

df <- df[,c('manuscript_id','Diagnosis','IPSS.M', 'Age.cat','IMPACT', 'SF3B1')]
colnames(df) <-  c('Study_ID','Diagnosis', 'IPSS-M','Age','IMPACT', 'SF3B1')
df$IMPACT[which(df$IMPACT == 'Not performed')] <- 'No'
df <- df[,rev(c('Study_ID','Diagnosis', 'IPSS-M','Age','IMPACT','SF3B1'))]
df[,colnames(df)] <- lapply(df[,colnames(df)] , factor)
df <- melt(df, id.var = c('Study_ID')) 

df$Study_ID <- factor(df$Study_ID, levels = c( "C-1", "C-2", "C-3", "C-4", "C-5", "C-6", "C-7", "C-8", "C-9",
                                               "C-10", "C-11", "C-12", "C-13", "C-14", "C-15", "C-16", "C-17", "C-18", "C-19", "C-20", "pA-1", "pA-2","pA-3","pA-4","pA-5",
  "MDS-1", "MDS-2", "MDS-3", "MDS-4", "MDS-5", "MDS-6", "MDS-7", "MDS-8", "MDS-9", "MDS-10",
  "MDS-11", "MDS-12", "MDS-13", "MDS-14", "MDS-15", "MDS-16", "MDS-17", "MDS-18", "MDS-19", "MDS-20",
  "MDS-21", "MDS-22", "MDS-23", "MDS-24", "MDS-25", "MDS-26", "MDS-27", "MDS-28", "MDS-29", "MDS-30",
  "MDS-31", "MDS-32", "MDS-33", "MDS-34", "MDS-35", "MDS-36", "MDS-37", "MDS-38", "MDS-39", "MDS-40",
  "MDS-41"))

df$value <- factor(df$value, levels = c('Control-Donor','Control-Lymph','Control-TP','MDS-LB','MDS-IB1','MDS-IB2', 'Very Low','Low','Moderate Low',
                                        'Moderate High', 'High', 'Very High', 'NA', 'Under 40','40-60', 'Over 60', 'Yes','No'))
brewer.pal(n = 8, name = "YlOrRd")


ggplot(df, aes(x = Study_ID, y = variable, fill = value)) +
  geom_tile(color = "black", lwd = 0.5) + 
  scale_fill_manual(values = c("No" = '#b5b5b5', "Yes" = "#405070", 
                               'NA' = "#b5b5b5",
                               'Very Low' = "#FFFFB2", 'Low' = "#FED976", 'Moderate Low' = "#FEB24C",
                               'Moderate High' = "#FD8D3C", 'High' ="#F03B20", 'Very High' = "#BD0026",
                               'Control-Donor' = "#0072B2", 'Control-Lymph' = "#6495ED",'Control-TP' = "#8CBADA",'MDS-LB' = "#ffd2c4", 'MDS-IB1' = "#FC9272", 'MDS-IB2' = "#DE2D26",
                               'Under 40'= '#c7d0ff','40-60'=  '#8d9ffc','Over 60' = '#3d4ea8')) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_grid(~ Study_ID, drop = TRUE, scales = "free", space = "free") + 
  theme_classic() + guides(fill=guide_legend(ncol=2)) +
  theme(
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    panel.spacing = unit(.05, "lines"),
    panel.border = element_blank(), 
    strip.background = element_blank(),
    strip.text = element_text(angle = 90, hjust = 0)) 


# --------------------------------------------------------------------------------
# Figure 1 panel c -- Cell-type frequency within samples, Control vs MDS (Xenium)
# --------------------------------------------------------------------------------

library(ggbreak)

df <- read.delim('data/Fig1c.csv', sep = ',')
df <- reshape2::melt(df)
df$diagnosis[which(df$diagnosis %in% c('Control-Donor','Control-Lymph'))] <- 'Control'
df$diagnosis[which(df$diagnosis %in% c('MDS-LB','MDS-IB1','MDS-IB2'))] <- 'MDS'

ggbarplot(
  df[which(df$diagnosis %in% c('Control','MDS')),], x = "variable", y = "value", 
  add = c("mean_sd", "jitter"), 
  color = "diagnosis", fill = 'diagnosis', palette = c("#0072B2","#BF504D"), alpha = 0.3,
  position = position_dodge(0.8), add.params = list(size = 0.5, alpha = 0.5),
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  stat_compare_means(aes(group = diagnosis), label = "p.signif", hide.ns=T, paired=F) + 
  labs(x = 'Cell Type', y = 'Proportion of Sample') + 
  ylim(c(-0.05, 0.8)) +
  scale_y_break(c(0.10, 0.25), space = 0.01, scales = 0.5)

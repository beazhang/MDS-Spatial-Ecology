# ============================================================================
# Extended Data Figure 2 -- R panels
#
# Run source("_setup.R") first for the shared library() calls.
# ============================================================================


# -----------------------------------------------------------------------------------
# Extended Data Figure 2 panel e -- Mean KIT expression in HSPCs, Control vs MDS-IB2 
# -----------------------------------------------------------------------------------

df <- read.delim('data/EDFig2e.csv', sep = ',')

ggplot(df, aes(x=diagnosis, y=KIT, color = diagnosis,fill = diagnosis)) + geom_boxplot(outliers = FALSE, alpha = 0.3) + theme_classic() + geom_jitter() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + 
  ylab('KIT Expression') + xlab('Condition') + scale_color_manual(values = c("#0072B2", "#a70000"), name = "Condition") + 
  scale_fill_manual(values = c("#0072B2", "#a70000"), name = 'Condition') + stat_compare_means(aes(group = diagnosis),hide.ns=T, paired=F) + ylim(0,1.5)

# ============================================================================
# Extended Data Figure 12 -- R panels
#
# Run source("_setup.R") first for the shared library() calls.
# ============================================================================


# --------------------------------------------------------------------------------------------------------------------
# Extended Data Figure 12 panel i -- Cytotoxic gene signature score of T cells in responders vs nonresponders (to HMA)
# --------------------------------------------------------------------------------------------------------------------

df <- read.delim('data/EDFig12i.csv', sep = ',')

ggplot(df, aes(x=cat, y=mean, color = cat,fill = cat)) + geom_boxplot(outliers = FALSE, alpha = 0.3, width = 0.75) + theme_classic() + geom_jitter() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + 
  ylab('Cytotoxic Gene Signature Score') + xlab('Condition') + scale_color_manual(values = c("#807F7F", "#BF504D"), name = "Condition") + 
  scale_fill_manual(values = c("#807F7F", "#BF504D"), name = 'Condition') + stat_compare_means(aes(group = cat),hide.ns=F, paired=F) + ylim(-0.5,0.5)



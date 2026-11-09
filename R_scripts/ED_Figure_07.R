# ============================================================================
# Extended Data Figure 7 -- R panels
#
# Run source("_setup.R") first for the shared library() calls.
# ============================================================================


# ------------------------------------------------------------------------------------------
# Extended Data Figure 7 panel c -- Mean VSIR / LGALS9 / S100A9 expression, Control vs MDS
# ------------------------------------------------------------------------------------------

### Monocytes
df <- read.delim('data/EDFig7c_Monocyte.csv',sep = ',')

ggplot(df, aes(x = type, y = VSIR, color = type, fill = type)) + ylab("VSIR") +
  geom_boxplot(alpha = 0.5, width = 0.6) + geom_point(size = 2, aes(color = type, alpha = 0.3)) +
  theme_classic() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.title = element_text(size=14), legend.text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(name = "Sample", values= (c("#0072B2",  "#a70000"))) + 
  scale_fill_manual(name = "Sample", values= (c("#0072B2",  "#a70000"))) + 
  stat_compare_means(label = "p.format",hide.ns=F, paired=F,  method.args = list(var.equal = FALSE)) + ylim(0.5, 2.5)

### Macrophage 
df <- read.delim('data/EDFig7c_Macrophage.csv',sep = ',')
ggplot(df, aes(x = type, y = LGALS9, color = type, fill = type)) + ylab("LGALS9") +
  geom_boxplot(alpha = 0.5, width = 0.6) + geom_point(size = 2, aes(color = type, alpha = 0.3)) +
  theme_classic() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.title = element_text(size=14), legend.text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(name = "Sample", values= (c("#0072B2",  "#a70000"))) + 
  scale_fill_manual(name = "Sample", values= (c("#0072B2",  "#a70000"))) + 
  stat_compare_means(label = "p.format",hide.ns=F, paired=F,  method.args = list(var.equal = FALSE)) + ylim(0, 1.5)


ggplot(df, aes(x = type, y = S100A9, color = type, fill = type)) + ylab("S100A9") +
  geom_boxplot(alpha = 0.5, width = 0.6) + geom_point(size = 2, aes(color = type, alpha = 0.3)) +
  theme_classic() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.title = element_text(size=14), legend.text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(name = "Sample", values= (c("#0072B2",  "#a70000"))) + 
  scale_fill_manual(name = "Sample", values= (c("#0072B2",  "#a70000"))) + 
  stat_compare_means(label = "p.format",hide.ns=F, paired=F, method.args = list(var.equal = FALSE)) + ylim(0, 2)


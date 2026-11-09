# ============================================================================
# Extended Data Figure 3 -- R panels
#
# Run source("_setup.R") first for the shared library() calls.
# ============================================================================


# ------------------------------------------------------------------------------------
# Extended Data Figure 3 panel b -- Cell-type proportion of sample by diagnosis group
# ------------------------------------------------------------------------------------

df <- read.delim('data/EDFig3b.csv', sep = ',')
df <- melt(df)
df$diagnosis[which(df$diagnosis %in% c('Control-Donor','Control-Lymph'))] <- 'Control'
df$diagnosis <- factor(df$diagnosis, levels = c('Control','MDS-LB','MDS-IB1','MDS-IB2'))
ggbarplot(
  df[which(df$diagnosis != 'Control-TP'),], x = "variable", y = "value", 
  add = c("mean_sd", "jitter"), 
  color = "diagnosis", fill = 'diagnosis', palette = c("#0072B2", "#FC9272", "#ed0e0e", "#a70000"), alpha = 0.3,
  position = position_dodge(0.8), add.params = list(size = 0.5, alpha = 0.5),
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Cell Type', y = 'Proportion of Sample') + ylim(c(-0.1, 0.8))


# ------------------------------------------------------------------------------------
# Extended Data Figure 3 panel c -- Lymphocyte subset frequency, Control vs MDS (Xenium)
# ------------------------------------------------------------------------------------

meta <- read.delim('data/EDFig3c.txt.gz', sep = '\t')
tmp <- meta[which(meta$ct_general %in% c('B Cell', 'T Cell','NK', 'pDC')),]
tmp <- as.data.frame(prop.table(table(tmp$ct_refined, tmp$study_id), 2))
tmp$type <- unlist(lapply(as.character(tmp$Var2), function(x) unlist(strsplit(x, '-'))[1]))
tmp$type[tmp$type %in% c('C')] <- 'Control'
tmp$type[tmp$type %in% c('MDS')] <- 'MDS'

tmp$Var1 <- factor(tmp$Var1, levels = c('B Immature', 'B Mature', 'CD4 Naive/Central Memory','CD4 Memory','CD4 Treg', 'CD8 Naive','CD8 Memory','NK', 'pDC'))
ggbarplot(
  tmp[tmp$type != 'pA',], x = "Var1", y = "Freq", 
  add = c("mean_sd", "jitter"),
  color = "type", fill = 'type', palette = c("#0072B2","#BF504D"), alpha = 0.3,
  position = position_dodge(0.8)
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + stat_compare_means(aes(group = type),label = "p.signif",hide.ns=T, paired=F) +
  ylab('Frequency within Lymphocytes') + xlab('Cell Type') + scale_color_manual(values = c("#0072B2","#BF504D"), name = "Type") +
  scale_fill_manual(values = c("#0072B2","#BF504D"), name = "Type") + scale_y_continuous(limits = c(-0, 0.5), breaks = seq(-0, 0.5, by = 0.1))

res <- compare_means(Freq ~ type, data = tmp[tmp$type != 'pA',], method = "wilcox.test", group.by = "Var1", p.adjust.method = 'BH')
res$p.signif.adj <- symnum(res$p.adj,corr = FALSE,na = FALSE,cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))

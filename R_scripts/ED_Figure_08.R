# ============================================================================
# Extended Data Figure 8 -- R panels
#
# Run source("_setup.R") first for the shared library() calls.
# ============================================================================


# ----------------------------------------------------------------------------------------
# Extended Data Figure 8 panel c -- Lymphocyte subset frequency, Control vs MDS (CITE-seq)
# ----------------------------------------------------------------------------------------

df <- read.delim('data/EDFig8c.csv', sep = ',')
df <- df[-which(df$type == 'MDS-Aza/Ven'),]
df <- reshape2::melt(df)
df$type[which(df$type %in% c('MDS-EB1','MDS-EB2'))] <- 'MDS'
df$type <- factor(df$type, levels = c('Control','MDS','MDS-EB1','MDS-EB2','MDS-Aza','MDS-PostAllo'))
df$variable <-  gsub('\\.', ' ', df$variable)
df$variable <- factor(df$variable, levels = c('B Immature', 'B Mature', 'CD4 Naive Central Memory', 'CD4 Treg','CD4 Memory', 'CD8 Naive', 'CD8 MAIT','CD8 Memory','CD8 gdT','NK BCL11B Bright', 'NK CD56 Bright', 'NK CD56 Dim'))
df <- df[-which(df$patient %in% c('ND44_T2','MDS_2_T3','ND41_T2')),]

ggbarplot(
  df[which(df$type %in% c('Control','MDS')),], x = "variable", y = "value",
  add = c("mean_sd", 'jitter'),
  color = "type", fill = 'type',  alpha = 0.3, palette = c('#0072B2', "#BF504D","#e08b89", "#BD0026", '#850310', '#a398e3', "#fcdefc"), 
  position = position_dodge(0.8)
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + stat_compare_means(aes(group = type),label = "p.signif",hide.ns=T, paired=F) +
  ylab('Frequency') + xlab('Cell Type') + scale_color_manual(values = c('#0072B2', "#BF504D","#e08b89", "#BD0026", '#850310', '#a398e3', "#fcdefc"), name = "type") + scale_y_continuous(limits = c(-0.2, 1), breaks = seq(-0.2, 1, by = 0.2))

res <- compare_means(value ~ type, data = df[which(df$type %in% c('Control','MDS')),], method = "wilcox.test", group.by = "variable", p.adjust.method = 'BH', paired = F)
res$p.signif.adj <- symnum(res$p.adj,corr = FALSE,na = FALSE,cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))


# --------------------------------------------------------------------------------------------------
# Extended Data Figure 8 panel d -- Proportion of CD4 and of CD8 T cells by subset, Control vs MDS
# --------------------------------------------------------------------------------------------------

### CD4
df <- read.delim('data/EDFig8d_CD4.csv', sep = ',') 
df$type[which(df$type %in% c('MDS-EB1','MDS-EB2'))] <- 'MDS'
df <- df[which(df$type %in% c('Control','MDS')),]
df <- reshape2::melt(df)

df$variable <- factor(df$variable, levels = c('CD4.Naive.Central.Memory','CD4.Memory','CD4.FOXO1.','CD4.Treg', 'CD8.Naive','CD8.gdT','CD8.MAIT','CD8.Effector.Memory','CD8.Memory.TOX.', 'CD8.Memory.GZMK.', 'CD8.TEMRA.GZMB.'))
df <- df[-which(df$patient %in% c('ND44_T2','MDS_2_T3','ND41_T2')),]
ggbarplot(
  df, x = "variable", y = "value", 
  add = c("mean_sd", "jitter"), 
  color = "type", fill = 'type', palette = c("#0072B2","#BF504D",'#a398e3'), alpha = 0.3,
  position = position_dodge(0.8)
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Cell Type', y = 'Proportion of CD4 T Cells') + ylim(c(-.1, 1))

res <- compare_means(value ~ type, data = df[which(df$type %in% c('Control','MDS')),], method = "wilcox.test", group.by = "variable", p.adjust.method = 'BH', paired = F)
res$p.signif.adj <- symnum(res$p.adj,corr = FALSE,na = FALSE,cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))
res

### CD8
df <- read.delim('data/EDFig8d_CD8.csv', sep = ',') 
df$type[which(df$type %in% c('MDS-EB1','MDS-EB2'))] <- 'MDS'
df <- df[which(df$type %in% c('Control','MDS')),]
df <- reshape2::melt(df)

df$variable <- factor(df$variable, levels = c('CD4.Naive.Central.Memory','CD4.Memory','CD4.FOXO1.','CD4.Treg', 'CD8.Naive','CD8.gdT','CD8.MAIT','CD8.Effector.Memory','CD8.Memory.TOX.', 'CD8.Memory.GZMK.', 'CD8.TEMRA.GZMB.'))
df <- df[-which(df$patient %in% c('ND44_T2','MDS_2_T3','ND41_T2')),]
ggbarplot(
  df, x = "variable", y = "value", 
  add = c("mean_sd", "jitter"), 
  color = "type", fill = 'type', palette = c("#0072B2","#BF504D",'#a398e3'), alpha = 0.3,
  position = position_dodge(0.8)
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Cell Type', y = 'Proportion of CD8 T Cells') + ylim(c(-.1, 1))

res <- compare_means(value ~ type, data = df[which(df$type %in% c('Control','MDS')),], method = "wilcox.test", group.by = "variable", p.adjust.method = 'BH', paired = F)
res$p.signif.adj <- symnum(res$p.adj,corr = FALSE,na = FALSE,cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))
res





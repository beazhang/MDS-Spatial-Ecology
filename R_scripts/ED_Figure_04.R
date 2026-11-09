# ============================================================================
# Extended Data Figure 4 -- R panels
#
# Run source("_setup.R") first for the shared library() calls.
# ============================================================================


# ------------------------------------------------------------------------------------------------------------------------------------------
# Extended Data Figure 4 panel b -- Cytotoxic gene signature score of T cells in vs out niche '1 T cell, 1 MSC, 1 Mk, 2 Late/Mature Myeloid'
# ------------------------------------------------------------------------------------------------------------------------------------------

df <- read.delim('data/EDFig4b.csv', sep = ',')

ggplot(df, aes(x=in_niche, y=mean, color = in_niche,fill = in_niche)) + geom_boxplot(outliers = FALSE, alpha = 0.3) + theme_classic() + geom_jitter(alpha = 0.3) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + 
  ylab('Cytotoxic Gene Signature') + xlab('Condition') + scale_color_manual(values = c("#0072B2", "#807F7F", "#a70000"), name = "Condition") + 
  scale_fill_manual(values = c("#0072B2",  "#807F7F","#a70000"), name = 'Condition') + stat_compare_means(ref.group = 'No', hide.ns=F, paired=T, label = "p.format") + ylim(-0.5,0.5)


# -------------------------------------------------------------------------------------------------------
# Extended Data Figure 4 panel d -- Niche enrichment for '1 Macrophage, 4 Erythroid' niche by diagnosis
# -------------------------------------------------------------------------------------------------------

currdf <- read.delim('data/Fig2e.csv.gz', sep = ',')  # shared input, named for Figure 2e
library(tidyr)
library(dplyr)
df_wide <- currdf[,c('sample', 'ratio','X', 'diagnosis', 'sample_id')] %>%
  pivot_wider(
    names_from = X, values_from = ratio
  )
df_wide <- as.data.frame(df_wide)
df_wide <- df_wide[, c('sample','22', 'diagnosis', 'sample_id')] # Niche 22: 4 Erythroid, 1 Macrophage
colnames(df_wide) <- c('sample','ratio','diagnosis', 'sample_id')
df_wide[is.na(df_wide)] <- 0

tmpdf <- df_wide
tmpdf <- tmpdf %>%
  group_by(sample_id, diagnosis) %>%
  summarise(
    ratio = mean(ratio, na.rm = TRUE),
    .groups = "drop"
  )
tmpdf$diagnosis[which(tmpdf$diagnosis %in% c('Control-Donor','Control-Lymph'))] <- 'Control'
tmpdf$diagnosis[which(tmpdf$diagnosis %in% c('Control-TP'))] <- 'Post Allo-HCT'
tmpdf$diagnosis <- factor(tmpdf$diagnosis, levels = c('Control','Post Allo-HCT','MDS-LB','MDS-IB1','MDS-IB2'))
ggplot(tmpdf, aes(x=diagnosis, y=ratio, color = diagnosis,fill = diagnosis)) + geom_boxplot(width=0.7, outlier.shape = NA, alpha = 0.3, position = position_dodge(0.9)) + theme_classic() + geom_jitter() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + 
  ylab('Niche Enrichment') + xlab('Condition') + scale_color_manual(values = c("#0072B2",  "#c4c4c4", "#FC9272", "#ed0e0e", "#a70000")) + scale_fill_manual(values = c("#0072B2",  "#c4c4c4", "#FC9272", "#ed0e0e", "#a70000")) 

wilcox.test(tmpdf$ratio[which(tmpdf$diagnosis == 'Control')], tmpdf$ratio[which(tmpdf$diagnosis %in% c('MDS-LB','MDS-IB1','MDS-IB2'))])


# --------------------------------------------------------------------------
# Extended Data Figure 4 panel g -- Immune aggregate size by diagnosis
# --------------------------------------------------------------------------

df <- read.delim('data/EDFig3c.txt.gz')  # shared input, named for Extended Data Figure 3c
df$type <- 'Control'
df$type[which(df$diagnosis %in% c('Control-TP'))] <- 'Post Allo-HCT'
df$type[which(df$diagnosis %in% c('MDS-IB1', 'MDS-IB2', 'MDS-LB'))] <- 'MDS'

tmp <- df[df$in_tls == 'Yes',]
tmp <- as.data.frame(table(tmp$tls,tmp$type))
tmp <- tmp[tmp$Freq > 0,]
tmp$Var2 <- factor(tmp$Var2, levels = c('Control','Post Allo-HCT','MDS'))
ggplot(tmp, aes(x=Var2, y=log10(Freq), color = Var2,fill = Var2)) + geom_boxplot(outliers = FALSE, alpha = 0.3) + theme_classic() + geom_jitter() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + 
  ylab('log10(Aggregate Size)') + xlab('Condition') + scale_color_manual(values = c("#0072B2", "#c4c4c4","#BF504D"), name = "Condition") + ggtitle("Immune Aggregate Size\n(# T, NK, B, pDC)") + 
  scale_fill_manual(values = c("#0072B2", "#c4c4c4","#BF504D"), name = 'Condition') + stat_compare_means(aes(group = Var2),label = "p.signif",hide.ns=T, paired=F)


# --------------------------------------------------------------------------
# Extended Data Figure 4 panel h -- Immune aggregate size by sample
# --------------------------------------------------------------------------

df <- read.delim('data/EDFig3c.txt.gz')  # shared input, named for Extended Data Figure 3c
df$type <- 'Control'
df$type[which(df$diagnosis %in% c('Control-TP'))] <- 'Post Allo-HCT'
df$type[which(df$diagnosis %in% c('MDS-IB1', 'MDS-IB2', 'MDS-LB'))] <- 'MDS'

tmp <- df[df$in_tls == 'Yes',]
tmp <- as.data.frame(table(tmp$tls,tmp$study_id))
tmp <- tmp[tmp$Freq > 0,]
tmp$Var2 <- factor(tmp$Var2, levels = c('C-1', 'C-2', 'C-3', 'C-4', 'C-5', 'C-6', 'C-7', 'C-8', 'C-9', 'C-10', 'C-11', 'C-12', 'C-13', 'C-14', 'C-15', 'pA-1', 'pA-2', 'pA-3', 'pA-4', 'pA-5', 'MDS-1', 'MDS-2', 'MDS-3', 'MDS-4', 'MDS-5', 'MDS-6', 'MDS-7', 'MDS-8', 'MDS-9', 'MDS-10', 'MDS-11', 'MDS-12', 'MDS-13', 'MDS-14', 'MDS-15', 'MDS-16', 'MDS-17', 'MDS-18', 'MDS-19', 'MDS-20', 'MDS-21', 'MDS-22', 'MDS-23', 'MDS-24', 'MDS-25', 'MDS-26', 'MDS-27', 'MDS-28', 'MDS-29', 'MDS-30', 'MDS-31', 'MDS-32', 'MDS-33', 'MDS-34', 'MDS-35', 'MDS-36', 'MDS-37', 'MDS-38', 'MDS-39', 'MDS-40', 'MDS-41'))
ggplot(tmp, aes(x=Var2, y=log10(Freq))) + geom_boxplot(outliers = FALSE, alpha = 0.3) + theme_classic() + geom_jitter(size = 0.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + 
  ylab('log10(Aggregate Size)') + xlab('Condition') + scale_color_manual(values = c("#0072B2", "#c4c4c4","#BF504D"), name = "Condition") + ggtitle("Immune Aggregate Size\n(# T, NK, B, pDC)") + 
  scale_fill_manual(values = c("#0072B2", "#c4c4c4","#BF504D"), name = 'Condition') 


# --------------------------------------------------------------------------
# Extended Data Figure 4 panel i -- Dominant cell type per immune aggregate
# --------------------------------------------------------------------------

df <- read.delim('data/EDFig3c.txt.gz')  # shared input, named for Extended Data Figure 3c
df <- df[which(df$in_tls == 'Yes'),]
df <- df[which(df$diagnosis != 'Control-TP'),]
df$type <- 'Control'
df$type[which(df$diagnosis %in% c('MDS-LB', 'MDS-IB1','MDS-IB2'))] <- 'MDS'
df_size <- as.data.frame(table(df$tls))
df_size$type <- df$type[match(df_size$Var1,df$tls)]
df_size$cat <- 'MDS Small'
df_size$cat[which(df_size$type == 'MDS' & df_size$Freq >= 100)] <- 'MDS Large'
df_size$cat[which(df_size$type == 'Control')] <- 'Control'

tmp <- as.data.frame(table(df$tls,df$ct_refined))
df_max <- tmp %>%
  group_by(Var1) %>%
  filter(Freq == max(Freq)) %>%
  ungroup()

df_max$type <- df_size$cat[match(df_max$Var1, df_size$Var1)]
df_plot <- melt(table(df_max$Var2, df_max$type))
df_plot$Var1 <- factor(df_plot$Var1, levels = c('B Immature', 'B Mature', 'CD4 Naive/Central Memory','CD4 Memory','CD4 Treg','CD8 Naive','CD8 Memory','NK', 'pDC'))
df_plot$Var2 <- factor(df_plot$Var2, levels = c('Control','MDS Small','MDS Large'))

df_plot <- df_plot[which(df_plot$value > 0),]
ggplot(df_plot, aes(fill=Var1, y=value, x=Var2)) + 
  geom_bar(position="fill", stat="identity",colour="black") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12), legend.text=element_text(size=12), legend.title=element_text(size=12)) +
  labs(x = 'Condition', y = 'Proportion of Immune Aggregates') +
  scale_fill_manual(name = "Dominating Cell Type", values= c("#F0E442", "#E69F00", "#56B4E9", "#0072B2", "#9cdba3", "#009E73","#f7b6d2"))


# --------------------------------------------------------------------------
# Extended Data Figure 4 panel k -- Lymphocyte frequency in immune aggregates
# --------------------------------------------------------------------------

df <- read.delim('data/EDFig3c.txt.gz')  # shared input, named for Extended Data Figure 3c
df <- df[which(df$in_tls == 'Yes'),]
df <- df[which(df$diagnosis != 'Control-TP'),]
df$type <- 'Control'
df$type[which(df$diagnosis %in% c('MDS-LB', 'MDS-IB1','MDS-IB2'))] <- 'MDS'
df_size <- as.data.frame(table(df$tls))
df_size$type <- df$type[match(df_size$Var1,df$tls)]
df_size$cat <- 'MDS Small'
df_size$cat[which(df_size$type == 'MDS' & df_size$Freq >= 100)] <- 'MDS Large'
df_size$cat[which(df_size$type == 'Control')] <- 'Control'

df_plot <- as.data.frame(prop.table(table(df$tls,df$ct_refined), margin = 1))
df_plot$type <- df_size$cat[match(df_plot$Var1, df_size$Var1)]

df_plot$Var2 <- factor(df_plot$Var2, levels = c('B Immature', 'B Mature', 'CD4 Naive/Central Memory','CD4 Memory','CD4 Treg','CD8 Naive','CD8 Memory','NK', 'pDC'))
df_plot$type <- factor(df_plot$type, levels = c('Control','MDS Small','MDS Large'))

ggbarplot(
  df_plot, x = "Var2", y = "Freq", #size = 0.001,
  add = c("mean_sd",'jitter'),
  color = "type", fill = 'type', palette = c("#0072B2","#e08b89","#BF504D"), alpha = 0.3,
  position = position_dodge(0.8), add.params = list(size = 0.5)
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  ylab('Frequency in Immune Aggregate') + xlab('Cell Type') + scale_color_manual(values = c("#0072B2","#e08b89", "#BF504D"), name = "type") 

res <- compare_means(Freq ~ type, data = df_plot, method = "wilcox.test", group.by = "Var2", p.adjust.method = 'BH')
res$p.signif.adj <- symnum(res$p.adj,corr = FALSE,na = FALSE,cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))



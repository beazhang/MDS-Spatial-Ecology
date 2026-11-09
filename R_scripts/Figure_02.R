# ============================================================================
# Figure 2 -- R panels
#
# Run source("_setup.R") first for the shared library() calls.
# ============================================================================


# ----------------------------------------------------------------------------------------------------------
# Figure 2 panel e -- Niche enrichment by diagnosis for niche '1 T cell, 1 MSC, 1 Mk, 2 Late/Mature Myeloid'
# ----------------------------------------------------------------------------------------------------------

currdf <- read.delim('data/Fig2e.csv.gz', sep = ',')
library(tidyr)
library(dplyr)
df_wide <- currdf[,c('sample', 'ratio','X', 'diagnosis', 'sample_id')] %>%
  pivot_wider(
    names_from = X, values_from = ratio
  )
df_wide <- as.data.frame(df_wide)
df_wide <- df_wide[, c('sample','490', 'diagnosis', 'sample_id')]
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
# Figure 2 panel f -- p-values for neighborhood enrichment heatmap
# --------------------------------------------------------------------------

df.tmp <- read.delim('data/Fig2f.txt', sep = '\t')

df.tmp$group <- paste(df.tmp$from, df.tmp$to)
res <- compare_means(enrichment ~ type, data = df.tmp, method = "wilcox.test", group.by = "group", p.adjust.method = 'BH')
res$p.signif.adj <- symnum(res$p.adj,corr = FALSE,na = FALSE,cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))


# ----------------------------------------------------------------------------------------
# Figure 2 panel g -- Erythroid-Macrophage neighborhood enrichment z-score control vs MDS
# ----------------------------------------------------------------------------------------

df <- read.delim('data/Fig2g.txt', sep = '\t')
df$type <- factor(df$type, levels = c('Control','MDS'))
df <- df[which(df$to == 'Macrophage'),]
df <- df[which(df$from == 'Erythroid'),]

ggbarplot(
  df, 
  x = "type", y = "enrichment", 
  add = c("mean_sd", "jitter"),
  color = "type", fill = "type", palette = c("#0072B2","#BF504D"), alpha = 0.3,
  position = position_dodge(0.8), add.params = list(size = 1), size = 1
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + guides(color = FALSE, fill = FALSE) +
  stat_compare_means(aes(group = type, label = paste0("p = ", after_stat(p.format)))) + ggtitle('Erythroid-Macrophage') + ylab('Enrichment z-score') + xlab('Condition') + scale_y_continuous(limits = c(-50, 25), breaks = seq(-50, 25, by = 25))


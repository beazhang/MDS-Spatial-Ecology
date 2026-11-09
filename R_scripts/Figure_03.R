# ============================================================================
# Figure 3 -- R panels
#
# Run source("_setup.R") first for the shared library() calls.
# ============================================================================


# ----------------------------------------------------------------------------
# Figure 3 panel b -- % of total T/NK/B/pDC in immune aggregates by diagnosis
# ----------------------------------------------------------------------------

df <- read.delim('data/Fig3b.csv', sep = ',')
df$Group <- factor(df$Group, levels = c('Control','Post Allo-HCT','MDS'))
ggplot(df, aes(x=Group, y=Value, color = Group,fill = Group)) + geom_jitter() + geom_boxplot(outliers = FALSE, alpha = 0.3) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + 
  ylab('% of T, NK, B, pDC\nin Immuune Aggregate') + xlab('Condition') + scale_color_manual(values = c("#0072B2", "#c4c4c4","#BF504D"), name = "Condition") + ggtitle("") + 
  scale_fill_manual(values = c("#0072B2", "#c4c4c4","#BF504D"), name = 'Condition') + 
  stat_compare_means(ref.group = 'MDS', hide.ns = TRUE, label = "p.format") + ylim(NA, 0.3) 

ggbarplot(
  df,
  x = "Group", y = "Value", 
  add = c("mean_sd", "jitter"),
  color = "Group", fill = 'Group', palette = c("#0072B2", "#c4c4c4","#BF504D"), alpha = 0.3,
  position = position_dodge(0.8)
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + stat_compare_means(ref.group = 'MDS', hide.ns = TRUE, label = "p.format") +
  ylab('% of T, NK, B, pDC\nin Immuune Aggregate') + xlab('Cell Type') + scale_color_manual(values = c("#0072B2", "#c4c4c4","#BF504D"), name = "Group") + ylim(NA, 0.3)


# -------------------------------------------------------------------------------
## Figure 3 panel f -- Volcano: MSCs in MDS, inside vs outside immune aggregates
# -------------------------------------------------------------------------------
library('ggrepel')

df <- read.delim('data/Fig3f.csv', sep = ',')

df$diffacessible <- 'Not Significant'
df$diffacessible[df$log2FoldChange > 0.25 & df$padj < 0.05] <- paste0("Up in Aggregate")
df$diffacessible[df$log2FoldChange < -0.25 & df$padj < 0.05] <- paste0("Up outside of Aggregate")
df$DAlabel <- df$X 
df$DAlabel[df$diffacessible == 'Not Significant'] <- NA
df$DAlabel[abs(df$log2FoldChange) < 0.25] <- NA
df$DAlabel[which(!(df$DAlabel %in% c('CCL19','CXCL9','CXCL10','CCL5','CXCL13')))] <- NA

mycolors <- c("#b2d0ed","#edb2b2", "gray")
names(mycolors) <- c(paste0("Up outside of Aggregate"), paste0("Up in Aggregate"), "Not Significant")

labelcolor <- df$diffacessible
labelcolor[which(df$diffacessible == 'Up outside of Aggregate')] <- '#01227d'
labelcolor[which(df$diffacessible == 'Up in Aggregate')] <- '#750000'
labelcolor[which(is.na(df$DAlabel))] <- NA

mycolors1 <- c('#01227d', "#b2d0ed","#750000", "#edb2b2","gray")
names(mycolors1) <- c(paste0("Out Aggregate Circle"), paste0("Up outside of Aggregate"), paste0("Aggregate Circle"), paste0("Up in Aggregate"), "Not Significant")
df$circle <- df$diffacessible 
df$circle[which((df$DAlabel != 'NA') & (df$log2FoldChange > 0))] <- 'Aggregate Circle'
df$circle[which((df$DAlabel != 'NA') & (df$log2FoldChange < 0))] <- 'Out Aggregate Circle'

ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj),label =DAlabel)) + 
  geom_point(aes(color = circle, fill = diffacessible), size = 3.1, stroke = 0.5, shape = 21) + 
  theme_classic() + geom_text_repel(size = 4, max.overlaps = Inf, color = labelcolor) + 
  xlab('log2(FC)') + ylab('-log10(adjusted p-value)') + 
  scale_fill_manual(name = "DE Isoforms", values = mycolors) + 
  scale_color_manual(name = "DE Isoforms", values = mycolors1) + theme(plot.title = element_text(hjust = 0.5))  + labs(title=paste0('MSCs in MDS')) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)) + scale_x_continuous(limits = c(-2.5, 5), breaks = seq(-2.5, 5, by = 2.5))


# --------------------------------------------------------------------------------
# Figure 3 panel g -- Volcano: T cells in MDS, inside vs outside immune aggregates
# --------------------------------------------------------------------------------
library('ggrepel')

df <- read.delim('data/Fig3g.csv', sep = ',')

df$diffacessible <- 'Not Significant'
df$diffacessible[df$log2FoldChange > 0.25 & df$padj < 0.05] <- paste0("Up in Aggregate")
df$diffacessible[df$log2FoldChange < -0.25 & df$padj < 0.05] <- paste0("Up outside of Aggregate")

df$DAlabel <- df$X 
df$DAlabel[df$diffacessible == 'Not Significant'] <- NA
df$DAlabel[abs(df$log2FoldChange) < 0.25] <- NA
df$DAlabel[which(!(df$DAlabel %in% c('CD244','KLRC1','IL7R', 'SELL','ENTPD1','KLRB1', 'IL2RA','TBX21','TCF7','CXCR4', 'CD27','CXCL9','CCL19','CXCL10','CXCL13', 'LEF1','KLRG1','PRF1','KLRD1','GZMA', 'GNLY', 'CTLA4','NKG7','CCL5','KLRD1','GZMK', 'IRF8', 'CCR7', 'CD8A','CD4', 'FOXP3')))] <- NA

mycolors <- c("#b2d0ed","#edb2b2", "gray")
names(mycolors) <- c(paste0("Up outside of Aggregate"), paste0("Up in Aggregate"), "Not Significant")

labelcolor <- df$diffacessible
labelcolor[which(df$diffacessible == 'Up outside of Aggregate')] <- '#01227d'
labelcolor[which(df$diffacessible == 'Up in Aggregate')] <- '#750000'
labelcolor[which(is.na(df$DAlabel))] <- NA

mycolors1 <- c('#01227d', "#b2d0ed","#750000", "#edb2b2","gray")
names(mycolors1) <- c(paste0("Out Aggregate Circle"), paste0("Up outside of Aggregate"), paste0("Aggregate Circle"), paste0("Up in Aggregate"), "Not Significant")
df$circle <- df$diffacessible
df$circle[which((df$DAlabel != 'NA') & (df$log2FoldChange > 0))] <- 'Aggregate Circle'
df$circle[which((df$DAlabel != 'NA') & (df$log2FoldChange < 0))] <- 'Out Aggregate Circle'

ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj),label =DAlabel)) + 
  geom_point(
    data = subset(df, diffacessible == "Not Significant"),
    aes(color = circle, fill = diffacessible),
    size = 3.1, stroke = 0.5, shape = 21
  ) +
  geom_point(
    data = subset(df, diffacessible == "Up outside of Aggregate"),
    aes(color = circle, fill = diffacessible),
    size =  2.75, stroke = 0.5, shape = 21
  ) +
  geom_point(
    data = subset(df, diffacessible == "Up in Aggregate"),
    aes(color = circle, fill = diffacessible),
    size =  2.75, stroke = 0.5, shape = 21
  ) +
  theme_classic() + geom_text_repel(size = 3, color=labelcolor, max.overlaps = Inf) + 
  xlab('log2(FC)') + ylab('-log10(adjusted p-value)') + 
  scale_fill_manual(name = "DE Isoforms", values = mycolors) + 
  scale_color_manual(name = "DE Isoforms", values = mycolors1) +
  theme(plot.title = element_text(hjust = 0.5)) + labs(title=paste0('T Cells in MDS')) + 
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, by = 3)) + scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5))

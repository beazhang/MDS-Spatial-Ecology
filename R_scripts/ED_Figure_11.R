# ============================================================================
# Extended Data Figure 11 -- R panels
#
# Run source("_setup.R") first for the shared library() calls.
# ============================================================================


# --------------------------------------------------------------------------------------
# Extended Data Figure 11 panel b -- Cells with >=1 ORAI2 canonical isoform transcript
# --------------------------------------------------------------------------------------

library(ggbreak) 
df <- read.delim('data/EDFig11b.csv', sep = ',')
df <- melt(df)
colnames(df) <- c('celltype','mutation','value')
df$celltype <- factor(df$celltype, levels = c("Adipocyte", "B Cell", "Early Myeloid", "Endothelial", "Erythroid", "HSPC", "Late/Mature Myeloid", "Macrophage", "Mast", "Megakaryocyte", "Mesenchymal/Stromal", "NK", "Osteoblast", "Plasma Cell", "T Cell", "VSMC", "pDC"))
ggplot(df, aes(x=celltype, y=value, color = mutation, fill = mutation)) + 
  geom_bar(stat = "identity",position=position_dodge2(padding = 0.25))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Cell Type', y = '# Cells') + scale_y_continuous(limits = c(0, 7000), breaks = seq(0, 7000, by = 1000)) + ggtitle('Canonical ORAI2 Cells') + scale_color_manual(values = c("#CC79A7","#009E73",  "#118cd4", "#E69F00", "#9aa0a6")) + scale_fill_manual(values = c("#CC79A7","#009E73",  "#118cd4", "#E69F00", "#9aa0a6"))


# -----------------------------------------------------------------------------------
# Extended Data Figure 11 panel c -- Cells with >=1 ORAI2 aberrant isoform transcript
# -----------------------------------------------------------------------------------

df <- read.delim('data/EDFig11c.csv', sep = ',')
df <- melt(df)
colnames(df) <- c('celltype','mutation','value')
df$celltype <- factor(df$celltype, levels = c("Adipocyte", "B Cell", "Early Myeloid", "Endothelial", "Erythroid", "HSPC", "Late/Mature Myeloid", "Macrophage", "Mast", "Megakaryocyte", "Mesenchymal/Stromal", "NK", "Osteoblast", "Plasma Cell", "T Cell", "VSMC", "pDC"))
ggplot(df, aes(x=celltype, y=value, color = mutation, fill = mutation)) + 
  geom_bar(stat = "identity",position=position_dodge2(padding = 0.25))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Cell Type', y = '# Cells') + scale_y_continuous(limits = c(0, 5500), breaks = seq(0, 5500, by = 500)) + 
  ggtitle('ORAI2 Mutant Cells') + scale_color_manual(values = c("#CC79A7","#009E73",  "#118cd4", "#E69F00", "#9aa0a6")) + scale_fill_manual(values = c("#CC79A7","#009E73",  "#118cd4", "#E69F00", "#9aa0a6"))


# ---------------------------------------------------------------------------
# Extended Data Figure 11 panel d -- Cells with >=2 SF3B1 mutant transcripts
# ---------------------------------------------------------------------------

df <- read.delim('data/EDFig11d.csv', sep = ',')
df <- melt(df)
colnames(df) <- c('celltype','mutation','value')
df$celltype <- factor(df$celltype, levels = c("Adipocyte", "B Cell", "Early Myeloid", "Endothelial", "Erythroid", "HSPC", "Late/Mature Myeloid", "Macrophage", "Mast", "Megakaryocyte", "Mesenchymal/Stromal", "NK", "Osteoblast", "Plasma Cell", "T Cell", "VSMC", "pDC"))
ggplot(df, aes(x=celltype, y=value, color = mutation, fill = mutation)) + 
  geom_bar(stat = "identity",position=position_dodge2(padding = 0.25))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Cell Type', y = '# Cells') + scale_y_continuous(limits = c(0, 1750), breaks = seq(0, 1750, by = 250)) + 
  ggtitle('>=2 SF3B1 Mutant Cells') + scale_color_manual(values = c("#CC79A7","#009E73",  "#118cd4", "#E69F00", "#9aa0a6")) + scale_fill_manual(values = c("#CC79A7","#009E73",  "#118cd4", "#E69F00", "#9aa0a6"))


# ---------------------------------------------------------------------------------------------------
# Extended Data Figure 11 panel e -- Cells with >=1 SF3B1 mutant and >=1 ZDHHC16 aberrant transcript
# ---------------------------------------------------------------------------------------------------

df <- read.delim('data/EDFig11e.csv', sep = ',')
df <- melt(df)
colnames(df) <- c('celltype','mutation','value')
df$celltype <- factor(df$celltype, levels = c("Adipocyte", "B Cell", "Early Myeloid", "Endothelial", "Erythroid", "HSPC", "Late/Mature Myeloid", "Macrophage", "Mast", "Megakaryocyte", "Mesenchymal/Stromal", "NK", "Osteoblast", "Plasma Cell", "T Cell", "VSMC", "pDC"))
ggplot(df, aes(x=celltype, y=value, color = mutation, fill = mutation)) + 
  geom_bar(stat = "identity",position=position_dodge2(padding = 0.25))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Cell Type', y = '# Cells') + scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10)) + 
  ggtitle('SF3B1 + ZDHHC16 Mutant Cells') + scale_color_manual(values = c("#CC79A7","#009E73",  "#118cd4", "#E69F00", "#9aa0a6")) + scale_fill_manual(values = c("#CC79A7","#009E73",  "#118cd4", "#E69F00", "#9aa0a6"))


# -----------------------------------------------------------------------------
# Extended Data Figure 11 panel g -- Cells with >=2 U2AF1 mutant transcripts
# -----------------------------------------------------------------------------

df <- read.delim('data/EDFig11g.csv', sep = ',')
df <- melt(df)
colnames(df) <- c('celltype','mutation','value')
df$celltype <- factor(df$celltype, levels = c("Adipocyte", "B Cell", "Early Myeloid", "Endothelial", "Erythroid", "HSPC", "Late/Mature Myeloid", "Macrophage", "Mast", "Megakaryocyte", "Mesenchymal/Stromal", "NK", "Osteoblast", "Plasma Cell", "T Cell", "VSMC", "pDC"))
ggplot(df, aes(x=celltype, y=value, color = mutation, fill = mutation)) + 
  geom_bar(stat = "identity",position=position_dodge2(padding = 0.25))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Cell Type', y = '# Cells') + scale_y_break(c(30, 70)) + scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, by = 10)) + 
  ggtitle('U2AF1 Mutant Cells') + scale_color_manual(values = c("#CC79A7","#009E73",  "#118cd4", "#E69F00", "#9aa0a6")) + scale_fill_manual(values = c("#CC79A7","#009E73",  "#118cd4", "#E69F00", "#9aa0a6"))


# ------------------------------------------------------------------------------------------------------------
# Extended Data Figure 11 panel h -- Proportion of U2AF1 mutant cells by cell type (in U2AF1 mutant samples)
# ------------------------------------------------------------------------------------------------------------

df <- read.delim('data/EDFig11h.csv', sep = ',')
df <- melt(df)
colnames(df) <- c('celltype','mutation','value')
ggplot(df[which(df$mutation == 'U2AF1'),], aes(x=celltype, y=value, color = mutation, fill = mutation)) + 
  geom_bar(stat = "identity",position=position_dodge2(padding = 1), width = 0.8)  + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Cell Type', y = '# Mutant Cells/# Total Cells') + scale_y_continuous(limits = c(0, 0.0015), breaks = seq(0, 0.0015, by = 0.0005)) + ggtitle('U2AF1 Mutant Cells') + 
  scale_color_manual(values = c( "#E69F00")) + scale_fill_manual(values = c("#E69F00"))


# ------------------------------------------------------------------------------------------------------------------------
# Extended Data Figure 11 panel j -- Frequency of T cell subtypes near vs far mutant/control myeloid/erythroid (alluvials)
# ------------------------------------------------------------------------------------------------------------------------

# Mutant Late/Mature Myeloid
df <- read.delim('data/EDFig11j_mutLateMyeloid.csv', sep = ',')

# Control Late/Mature Myeloid
#df <- read.delim('data/EDFig11j_ctlLateMyeloid.csv', sep = ',')

# Mutant Erythroid
#df <- read.delim('data/EDFig11j_mutErythroid.csv', sep = ',')

# Control Erythroid
#df <- read.delim('data/EDFig11j_ctlErythroid.csv', sep = ',')

df <- melt(df)
df$variable <- factor(df$variable, levels = rev(c("CD4.Naive.Central.Memory","CD4.Treg","CD4.Memory","CD8.Naive","CD8.Memory")))
ggplot(df, aes(x = near_mut, y = value, stratum = variable, alluvium = variable, fill = variable)) + ggalluvial::geom_flow() + theme_classic() +
  geom_stratum(alpha = .8) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values = (c("#009E73", '#d3f5a9', "#0072B2", "#68a2ed","#a0daf2"))) 

# ----------------------------------------------------------------------------------------------------------
# Extended Data Figure 11 panel k -- CD8 Memory T frequency (out of all T cells) near vs far from mutant Mk
# ----------------------------------------------------------------------------------------------------------

df <- read.delim('data/EDFig11k.csv', sep = ',')
df <- df[which(df$patient_id != 'ND170'),]
ggpaired(
  df, x = "near_mut", y = "CD8.Memory", alpha = 0.3,
  color = "near_mut", palette = c("#807F7F", "#BF504D"),alpha = 0.3,
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  stat_compare_means(aes(group = near_mut),hide.ns=T, paired=T) + ylab('Frequency') + xlab('Near Mutant') + 
  scale_color_manual(values = c("#807F7F", "#BF504D"), name = "Near Mutant") + scale_fill_manual(values = c("#807F7F", "#BF504D"), name = "Near Mutant") + 
scale_y_continuous(limits = c(0.0, 0.75), breaks = seq(0.0, 0.75, by = 0.25))


# ----------------------------------------------------------------------------
# Extended Data Figure 11 panel l -- Volcano: mutant vs control Megakaryocyte
# ----------------------------------------------------------------------------
library('ggrepel')

df <- read.delim('data/EDFig11l.csv', sep = ',')

df$diffacessible <- 'Not Significant'
df$diffacessible[df$log2FoldChange > 0.25 & df$padj < 0.05] <- paste0("Up in Mutant")
df$diffacessible[df$log2FoldChange < -0.25 & df$padj < 0.05] <- paste0("Up in Control")

df$DAlabel <- df$X 
df$DAlabel[df$diffacessible == 'Not Significant'] <- NA
df$DAlabel[which(!(df$DAlabel %in% c("TGFB1", "GP1BA","PECAM1","CCL5","CD4","IFITM3",
                                     'LAMP1','CXCR4','FCGR2A','SSR4','PSMB10','APOE')))] <- NA

mycolors <- c("#b2d0ed","#edb2b2", "gray")
names(mycolors) <- c(paste0("Up in Control"), paste0("Up in Mutant"), "Not Significant")

labelcolor <- df$diffacessible
labelcolor[which(df$diffacessible == 'Up in Control')] <- '#01227d'
labelcolor[which(df$diffacessible == 'Up in Mutant')] <- '#750000'
labelcolor[which(is.na(df$DAlabel))] <- "NA"

mycolors1 <- c('#01227d', "#b2d0ed","#750000", "#edb2b2","gray")
names(mycolors1) <- c(paste0("Down Circle"), paste0("Up in Control"), paste0("Up Circle"), paste0("Up in Mutant"), "Not Significant")
df$circle <- df$diffacessible
df$circle[which((df$DAlabel != 'NA') & (df$log2FoldChange > 0))] <- 'Up Circle'
df$circle[which((df$DAlabel != 'NA') & (df$log2FoldChange < 0))] <- 'Down Circle'

ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj),label =DAlabel)) + 
  geom_point(
    data = subset(df, diffacessible == "Not Significant"),
    aes(color = circle, fill = diffacessible),
    size = 3.1, stroke = 0.5, shape = 21
  ) +
  geom_point(
    data = subset(df, diffacessible == "Down near MSC"),
    aes(color = circle, fill = diffacessible),
    size = 3.1, stroke = 0.5, shape = 21
  ) +
  geom_point(
    data = subset(df, diffacessible == "Up in Mutant" & is.na(DAlabel)),
    aes(color = circle, fill = diffacessible),
    size = 3.1, stroke = 0.5, shape = 21
  ) +
  geom_point(
    data = subset(df, diffacessible == "Up in Mutant" & !is.na(DAlabel)),
    aes(color = circle, fill = diffacessible),
    size = 3.1, stroke = 0.5, shape = 21
  ) +
  theme_classic() + geom_text_repel(size = 3, color=labelcolor, max.overlaps = Inf) + 
  xlab('log2(FC)') + ylab('-log10(adjusted p-value)') + 
  scale_fill_manual(name = "DE Isoforms", values = mycolors) + 
  scale_color_manual(name = "DE Isoforms", values = mycolors1) +
  theme(plot.title = element_text(hjust = 0.5)) + labs(title=paste0('Mutant Mk vs Control Mk')) + 
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1)) + scale_x_continuous(limits = c(-4, 3), breaks = seq(-4, 3, by = 1))


# --------------------------------------------------------------------------
# Extended Data Figure 11 panel m -- Volcano: Megakaryocyte near vs far MSC
# --------------------------------------------------------------------------
library('ggrepel')

df <- read.delim('data/EDFig11m.csv', sep = ',')

df$diffacessible <- 'Not Significant'
df$diffacessible[df$log2FoldChange > 0.25 & df$padj < 0.05] <- paste0("Up near MSC")
df$diffacessible[df$log2FoldChange < -0.25 & df$padj < 0.05] <- paste0("Down near MSC")

df$DAlabel <- df$X 
df$DAlabel[df$diffacessible == 'Not Significant'] <- NA
df$DAlabel[which(!(df$DAlabel %in% c("TGFB1", "ACTA2", "IGFBP7", "CDKN2D", "CCND1",
                                     "JAK1", "JAK2", "JAK3", "STAT5B", "STAT4",
                                     "NOTCH1", "NOTCH3", "HIF1A", "SOCS2", "TNF",
                                     "CD44", "ITGA1", "IRF1", 'CXCL12')))] <- NA

mycolors <- c("#b2d0ed","#edb2b2", "gray")
names(mycolors) <- c(paste0("Down near MSC"), paste0("Up near MSC"), "Not Significant")

labelcolor <- df$diffacessible
labelcolor[which(df$diffacessible == 'Down near MSC')] <- '#01227d'
labelcolor[which(df$diffacessible == 'Up near MSC')] <- '#750000'
labelcolor[which(is.na(df$DAlabel))] <- "NA"

mycolors1 <- c('#01227d', "#b2d0ed","#750000", "#edb2b2","gray")
names(mycolors1) <- c(paste0("Down Circle"), paste0("Down near MSC"), paste0("Up Circle"), paste0("Up near MSC"), "Not Significant")
df$circle <- df$diffacessible
df$circle[which((df$DAlabel != 'NA') & (df$log2FoldChange > 0))] <- 'Up Circle'
df$circle[which((df$DAlabel != 'NA') & (df$log2FoldChange < 0))] <- 'Down Circle'

ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj),label =DAlabel)) + 
  geom_point(
    data = subset(df, diffacessible == "Not Significant"),
    aes(color = circle, fill = diffacessible),
    size = 3.1, stroke = 0.5, shape = 21
  ) +
  geom_point(
    data = subset(df, diffacessible == "Down near MSC"),
    aes(color = circle, fill = diffacessible),
    size = 3.1, stroke = 0.5, shape = 21
  ) +
  geom_point(
    data = subset(df, diffacessible == "Up near MSC" & is.na(DAlabel)),
    aes(color = circle, fill = diffacessible),
    size = 3.1, stroke = 0.5, shape = 21
  ) +
  geom_point(
    data = subset(df, diffacessible == "Up near MSC" & !is.na(DAlabel)),
    aes(color = circle, fill = diffacessible),
    size = 3.1, stroke = 0.5, shape = 21
  ) +
  theme_classic() + geom_text_repel(size = 3, color=labelcolor, max.overlaps = Inf) + 
  xlab('log2(FC)') + ylab('-log10(adjusted p-value)') + 
  scale_fill_manual(name = "DE Isoforms", values = mycolors) + 
  scale_color_manual(name = "DE Isoforms", values = mycolors1) +
  theme(plot.title = element_text(hjust = 0.5)) + labs(title=paste0('Mk near vs far MSC')) + 
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 5)) + scale_x_continuous(limits = c(-3, 2.), breaks = seq(-3, 2., by = 1))


# ----------------------------------------------------------------------------------
# Extended Data Figure 11 panel n -- Volcano: mutant vs control Late/Mature Myeloid
# ----------------------------------------------------------------------------------
library('ggrepel')

df <- read.delim('data/EDFig11n.csv', sep = ',')

df$diffacessible <- 'Not Significant'
df$diffacessible[df$log2FoldChange > 0.25 & df$padj < 0.05] <- paste0("Up in Mutant")
df$diffacessible[df$log2FoldChange < -0.25 & df$padj < 0.05] <- paste0("Up in Control")

df$DAlabel <- df$X 
df$DAlabel[df$diffacessible == 'Not Significant'] <- NA
df$DAlabel[abs(df$log2FoldChange) < 0.25] <- NA
df$DAlabel[which(!(df$DAlabel %in% c('ARG1','S100A9','CEACAM1')))] <- NA

mycolors <- c("#b2d0ed","#edb2b2", "gray")
names(mycolors) <- c(paste0("Up in Control"), paste0("Up in Mutant"), "Not Significant")

labelcolor <- df$diffacessible
labelcolor[which(df$diffacessible == 'Up in Control')] <- '#01227d'
labelcolor[which(df$diffacessible == 'Up in Mutant')] <- '#750000'
labelcolor[which(is.na(df$DAlabel))] <- NA

mycolors1 <- c('#01227d', "#b2d0ed","#750000", "#edb2b2","gray")
names(mycolors1) <- c(paste0("Down Circle"), paste0("Up in Control"), paste0("Up Circle"), paste0("Up in Mutant"), "Not Significant")
df$circle <- df$diffacessible
df$circle[which((df$DAlabel != 'NA') & (df$log2FoldChange > 0))] <- 'Up Circle'
df$circle[which((df$DAlabel != 'NA') & (df$log2FoldChange < 0))] <- 'Down Circle'


ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj),label =DAlabel)) + 
  geom_point(
    data = subset(df, diffacessible == "Not Significant"),
    aes(color = circle, fill = diffacessible),
    size = 3.1, stroke = 0.5, shape = 21
  ) +
  geom_point(
    data = subset(df, diffacessible == "Up in Control"),
    aes(color = circle, fill = diffacessible),
    size = 3.1, stroke = 0.5, shape = 21
  ) +
  geom_point(
    data = subset(df, diffacessible == "Up in Mutant" & is.na(DAlabel)),
    aes(color = circle, fill = diffacessible),
    size = 3.1, stroke = 0.5, shape = 21
  ) +
  geom_point(
    data = subset(df, diffacessible == "Up in Mutant" & !is.na(DAlabel)),
    aes(color = circle, fill = diffacessible),
    size = 3.1, stroke = 0.5, shape = 21
  ) +
  theme_classic() + geom_text_repel(size = 3, color=labelcolor, max.overlaps = Inf) + 
  xlab('log2(FC)') + ylab('-log10(adjusted p-value)') + 
  scale_fill_manual(name = "DE Isoforms", values = mycolors) + 
  scale_color_manual(name = "DE Isoforms", values = mycolors1) +
  theme(plot.title = element_text(hjust = 0.5)) + labs(title=paste0('Mutant vs Control Late/Mature Myeloid')) + 
  scale_x_continuous(limits = c(-4,3), breaks = seq(-4,3, by = 1)) + scale_y_continuous(limits = c(0,10), breaks = seq(0, 10, by = 2))



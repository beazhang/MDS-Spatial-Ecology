# ============================================================================
# Figure 4 -- R panels
#
# Run source("_setup.R") first for the shared library() calls.
# ============================================================================


# --------------------------------------------------------------------------
# Figure 4 panel h -- SF3B1 mutant cell counts by cell type
# --------------------------------------------------------------------------
library('ggbreak') 

df <- read.delim('data/Fig4h.csv', sep = ',')
df <- melt(df)
colnames(df) <- c('celltype','mutation','value')
df$celltype <- factor(df$celltype, levels = c("Adipocyte", "B Cell", "Early Myeloid", "Endothelial", "Erythroid", "HSPC", "Late/Mature Myeloid", "Macrophage", "Mast", "Megakaryocyte", "Mesenchymal/Stromal", "NK", "Osteoblast", "Plasma Cell", "T Cell", "VSMC", "pDC"))
ggplot(df, aes(x=celltype, y=value, color = mutation, fill = mutation)) + 
  geom_bar(stat = "identity",position=position_dodge2(padding = 0.25))  + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Cell Type', y = '# Cells') + scale_y_break(c(3000, 5000))  + scale_y_continuous(limits = c(0, 5500), breaks = seq(0, 5500, by = 500)) + ggtitle('SF3B1 Mutant Cells') +
  scale_color_manual(values = c("#CC79A7","#009E73",  "#118cd4", "#E69F00", '#9aa0a6')) + scale_fill_manual(values = c("#CC79A7","#009E73",  "#118cd4", "#E69F00",'#9aa0a6'))


# ---------------------------------------------------------------------------------------------
# Figure 4 panel i -- Proportion of SF3B1 mutant cells by cell type (in SF3B1 mutant samples)
# ---------------------------------------------------------------------------------------------

df <- read.delim('data/Fig4i.csv', sep = ',')
df <- melt(df)
colnames(df) <- c('celltype','mutation','value')
df$celltype <- factor(df$celltype, levels = c("Adipocyte", "B Cell", "Early Myeloid", "Endothelial", "Erythroid", "HSPC", "Late/Mature Myeloid", "Macrophage", "Mast", "Megakaryocyte", "Mesenchymal/Stromal", "NK", "Osteoblast", "Plasma Cell", "T Cell", "VSMC", "pDC"))
ggplot(df[which(df$mutation == 'SF3B1'),], aes(x=celltype, y=value, color = mutation, fill = mutation)) + 
  geom_bar(stat = "identity",position=position_dodge2(padding = 1), width = 0.8)  + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  labs(x = 'Cell Type', y = '# Mutant Cells/# Total Cells') + scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by = 0.025)) + ggtitle('SF3B1 Mutant Cells') + 
  scale_color_manual(values = c( "#118cd4")) + scale_fill_manual(values = c("#118cd4"))

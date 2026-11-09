# ============================================================================
# Extended Data Figure 10 -- R panels
#
# Run source("_setup.R") first for the shared library() calls.
# ============================================================================


# ---------------------------------------------------------------------------------------------
# Extended Data Figure 10 panel c -- TRAC and AHSP expression in erythroid cells by TCR status
# ---------------------------------------------------------------------------------------------

df <- read.delim('data/EDFig10c.txt.gz', sep = '\t')

ggplot(df, aes(x = TCRpos, y = AHSP, color = TCRpos, fill = TCRpos)) + ylab("AHSP Expression") +
  geom_violin(linewidth = 0.8, scale="width", alpha = 0.3,position = position_dodge(0.9)) + 
  theme_classic() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.title = element_text(size=14), legend.text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(name = "Sample", values=(c("#807F7F", "#BF504D"))) + 
  scale_fill_manual(name = "Sample", values= (c("#807F7F", "#BF504D"))) + 
  stat_compare_means(aes(group = TCRpos),hide.ns=F)  +scale_y_continuous(limits = c(0,4), breaks = seq(0,4, by = 1))


# ------------------------------------------------------------------------------------------------------------
# Extended Data Figure 10 panels g & l -- Cell-type composition of global vs TCR+ neighborhoods (MDS-2 T1, T3)
# ------------------------------------------------------------------------------------------------------------

## alluvials for each sample (example MDS-2 T4)
df <- read.delim('data/EDFig10g_l.txt', sep = ',')
df$ct_general <- factor(df$ct_general, levels = rev(c("Adipocyte", "B Cell", "Early Myeloid", "Endothelial", "Erythroid", "HSPC", "Late/Mature Myeloid", "Macrophage", "Mast", "Megakaryocyte", "Mesenchymal/Stromal", "NK", "Osteoblast", "Plasma Cell", "T Cell", "VSMC", "pDC")))

ggplot(df[df$patient_id == 'ND72_T4',], aes(x = type, y = value, stratum = ct_general, alluvium = ct_general, fill = ct_general)) + ggalluvial::geom_flow() + theme_classic() +
  geom_stratum(alpha = .8) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values = rev(c("#1f77b4", "#ff7f0e", "#4db073","#d62728", "#c5b0d5",'#8c564b', "#e377c2", "#035239", "#17becf", "#aec7e8", "#ffbb78",'#615db0',"#98df8a", "#ff9896","#a42aeb","#c49c94", "#f7b6d2"))) 


# ------------------------------------------------------------------------------------------------------------
# Extended Data Figure 10 panels i -- T cell subtype composition of cells with targeted TCRs (CITE vs Xenium)
# ------------------------------------------------------------------------------------------------------------

## TCR+ cell types (alluvial comparing CITE with xenium)
df <- read.delim('data/EDFig10i.txt', sep =',')
df$ct_refined[which(df$ct_refined == 'CD8 TEMRA GZMB+')] <- 'CD8 Memory'
df <- df %>%
  group_by(ct_refined) %>%
  summarise(
    xen = sum(xen, na.rm = TRUE),
    sc = sum(sc, na.rm = TRUE)
  )

df$xen <- df$xen/sum(df$xen)
df$sc <- df$sc/sum(df$sc)
df <- melt(df)
df$ct_refined <- factor(df$ct_refined, levels = c('CD8 Memory', 'CD8 TEMRA GZMB+',  'CD8 Naive/Central Memory', 'gdT','NK','CD4 Memory', 'CD4 Treg','CD4 Naive/Central Memory'))
df$variable <- factor(df$variable, levels = c('sc','xen'))

ggplot(df, aes(x = variable, y = value, stratum = ct_refined, alluvium = ct_refined, fill = ct_refined)) + ggalluvial::geom_flow() + theme_classic() +
  geom_stratum(alpha = .8) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values = c("#009E73", '#d3f5a9','#ffd13b',"#615db0", "#0072B2", "#68a2ed","#a0daf2"))


# -----------------------------------------------------------------------------------
# Extended Data Figure 10 panel j -- T cell frequency in global vs TCR+ neighborhood
# -----------------------------------------------------------------------------------

df <- read.delim('data/EDFig10j.txt', sep = ',')
df$type <- factor(df$type, levels = c('Global','TCR'))

ggpaired(
  df[which(df$ct_general == 'T Cell'),], 
  x = "type", y = "value", alpha = 0.3,
  color = "type", palette = c("#807F7F", "#BF504D"),alpha = 0.3,
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  stat_compare_means(aes(group = type),hide.ns=T, paired=T) + ylab('Frequency within Neighborhood') + xlab('Neighborhood') + 
  scale_color_manual(values = c("#807F7F", "#BF504D"), name = "Neighborhood") + scale_fill_manual(values = c("#807F7F", "#BF504D"), name = "Neighborhood") + 
  scale_y_continuous(limits = c(0.0, 0.075), breaks = seq(0.0, 0.075, by = 0.025))


# ------------------------------------------------------------------------------------------------
# Extended Data Figure 10 panel m -- CD8 memory frequency within immune aggregate, TCR- vs TCR+
# ------------------------------------------------------------------------------------------------

df <- read.delim('data/EDFig10m.csv', sep = ',')
df$tcr[which(df$tcr != 'None')] <- 'TCR+'
df$tcr[which(df$tcr == 'None')] <- 'TCR-'
df$tcr <- factor(df$tcr, levels = c('TCR-', 'TCR+'))
ggpaired(
  df, 
  x = "tcr", y = "Yes", alpha = 0.3,
  color = "tcr", palette = c("#807F7F", "#BF504D"),alpha = 0.3,
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  stat_compare_means(aes(group = tcr),hide.ns=T, paired=T) + ylab('Frequency within Immune Aggregate') + xlab('TCR') + 
  scale_color_manual(values = c("#807F7F", "#BF504D"), name = "TCR") + scale_fill_manual(values = c("#807F7F", "#BF504D"), name = "TCR") + scale_y_continuous(limits = c(0.0, 0.06), breaks = seq(0.0, 0.06, by = 0.02))




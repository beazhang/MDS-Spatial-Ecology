# ============================================================================
# Figure 5 -- R panels
#
# Run source("_setup.R") first for the shared library() calls.
# ============================================================================


# ------------------------------------------------------------------------------------------------------------------------
# Figure 5 panel a -- Frequency of T cell subtypes near vs far mutant/control Megakaryocytes (alluvials)
# ------------------------------------------------------------------------------------------------------------------------

# Mutant Megakaryocyte
df <- read.delim('data/Fig5a_mutMk.csv', sep = ',')

# Control Megakaryocyte
#df <- read.delim('data/Fig5a_ctlMk.csv', sep = ',')

df <- melt(df)
df$variable <- factor(df$variable, levels = rev(c("CD4.Naive.Central.Memory","CD4.Treg","CD4.Memory","CD8.Naive","CD8.Memory")))
ggplot(df, aes(x = near_mut, y = value, stratum = variable, alluvium = variable, fill = variable)) + ggalluvial::geom_flow() + theme_classic() +
  geom_stratum(alpha = .8) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values = (c("#009E73", '#d3f5a9', "#0072B2", "#68a2ed","#a0daf2")))


# ---------------------------------------------------------------------------------
# Figure 5 panel b -- Mean Mk TGFB1 expression, Control vs SF3B1 mutant (Xenium v2)
# ---------------------------------------------------------------------------------

df <- read.delim('data/Fig5b.csv',sep = ',')
ggplot(df, aes(x = mutant, y = TGFB1, color = mutant, fill = mutant)) + ylab("TGFB1") +
  geom_boxplot(alpha = 0.5, width = 0.75) + geom_point(size = 2, aes(color = mutant, alpha = 0.3)) +
  theme_classic() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.title = element_text(size=14), legend.text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(name = "Sample", values= (c("#807F7F", "#BF504D"))) + 
  scale_fill_manual(name = "Sample", values= (c("#807F7F", "#BF504D"))) + 
  stat_compare_means(label = "p.format",hide.ns=F, paired=F, method = 't.test', method.args = list(var.equal = FALSE)) + ylim(0, 3) 


# -----------------------------------------------------------------------
# Figure 5 panel c -- Mean Mk TGFB1 expression by diagnosis (Xenium v1)
# -----------------------------------------------------------------------

df <- read.delim('data/Fig5c.csv',sep = ',')
df$diagnosis[df$diagnosis %in% c('Control-Donor', 'Control-Lymph')] <- 'Control'
df$diagnosis <- factor(df$diagnosis, levels = c('Control', 'Control-TP', 'MDS-LB', 'MDS-IB1','MDS-IB2'))

ggplot(df[which(df$diagnosis != 'Control-TP'),], aes(x = diagnosis, y = TGFB1, color = diagnosis, fill = diagnosis)) + ylab("TGFB1") +
  geom_boxplot(alpha = 0.5, width = 0.5) + geom_point(size = 2, aes(color = diagnosis, alpha = 0.3)) +
  theme_classic() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14), legend.title = element_text(size=14), legend.text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(name = "Sample", values= (c("#0072B2", "#FC9272", "#ed0e0e", "#a70000"))) + 
  scale_fill_manual(name = "Sample", values= (c("#0072B2", "#FC9272", "#ed0e0e", "#a70000"))) + 
  stat_compare_means(label = "p.format",ref.group = 'MDS-IB2',hide.ns=F, paired=F) + ylim(0, 2.5)


# -------------------------------------------------------------------------------
# Figure 5 panel e -- Mean Mk TGFB1 expression near vs not near immune aggregate
# -------------------------------------------------------------------------------

df <- read.delim('data/Fig5e.csv',sep = ',')
keepidx <- c("MDS-11", "MDS-14", "MDS-15", "MDS-18", "MDS-20", "MDS-26", "MDS-28", "MDS-31", "MDS-35", "MDS-36", "MDS-41") # samples with adequate quantity of Mk and TGFB1 in xenium panel
df <- df[which(df$study_id %in% keepidx),]
df <- df[which(df$diagnosis %in% c('MDS-LB', 'MDS-IB1','MDS-IB2')), ]
ggpaired(df, x = "near_TLS", y = "TGFB1", alpha = 0.3, color = "near_TLS", group = 'sample_id', palette = c("#807F7F", "#BF504D"),
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  stat_compare_means(aes(group = near_TLS),hide.ns=T, paired=T) + ylab('TGFB1 Expression') + xlab('Near TLS') + 
  scale_color_manual(values = c("#807F7F", "#BF504D"), name = "Near TLS") + scale_fill_manual(values = c("#807F7F", "#BF504D"), name = "Near TLS") +
  scale_y_continuous(limits = c(0.5, 3.5), breaks = seq(0.5, 3.5, by = 0.5))


# ----------------------------------------------------------------------------------------------------
# Figure 5 panel h -- Cytotoxic gene signature score of CD8 Memory T cells by distance to nearest Mk
# ----------------------------------------------------------------------------------------------------

df <- read.delim('data/Fig5h.csv', sep = ',')

ggplot(df, aes(x=cat, y=mean, color = cat,fill = cat)) + geom_boxplot(outliers = FALSE, alpha = 0.3) + theme_classic() + geom_jitter(alpha = 0.3) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5), axis.text=element_text(size=12)) + 
  ylab('Cytotoxic Gene Signature') + xlab('Condition') + scale_color_manual(values = c("#0072B2", "#807F7F", "#a70000"), name = "Condition") + 
  scale_fill_manual(values = c("#0072B2",  "#807F7F","#a70000"), name = 'Condition') + stat_compare_means(ref.group = 'Far', hide.ns=T, paired=T, label = "p.format") + ylim(-0.6,0.6)



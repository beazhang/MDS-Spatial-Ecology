# ============================================================================
# Extended Data Figure 9 -- R panels
#
# Run source("_setup.R") first for the shared library() calls.
# ============================================================================


# ----------------------------------------------------------------------------------------------------------
# Extended Data Figure 9 panel a -- CD45dim subset frequency, Control vs MDS-PostAllo vs MDS-Aza (CITE-seq)
# ----------------------------------------------------------------------------------------------------------

df <- read.delim('data/EDFig9a.csv', sep = ',')
df <- df[-which(df$type == 'MDS-Aza/Ven'),]
df <- reshape2::melt(df)
df$type[which(df$type %in% c('MDS-EB1','MDS-EB2'))] <- 'MDS'
df$type <- factor(df$type, levels = c('Control','MDS','MDS-EB1','MDS-EB2','MDS-Aza','MDS-PostAllo'))
df$variable <-  gsub('\\.', ' ', df$variable)
df <- df[-which(df$patient %in% c('ND44_T2','MDS_2_T3','ND41_T2')),]
df$type <- factor(df$type, levels = c('Control','MDS-PostAllo','MDS','MDS-EB1','MDS-EB2', 'MDS-Aza'))
ggbarplot(
  df[which(df$type %in% c('Control','MDS-PostAllo', 'MDS-Aza')),], x = "variable", y = "value", #size = 0.001,
  add = c("mean_sd", 'jitter'),
  color = "type", fill = 'type',  alpha = 0.3, palette = c('#0072B2', "#fcdefc", '#a398e3'), 
  position = position_dodge(0.8)
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  ylab('Frequency') + xlab('Cell Type') + scale_color_manual(values = c('#0072B2', "#fcdefc", '#a398e3'), name = "type") + scale_y_continuous(limits = c(-0.2, 1), breaks = seq(-0.2, 1, by = 0.2))

res <- compare_means(value ~ type, data = df[which(df$type %in% c('Control','MDS-PostAllo', 'MDS-Aza')),], method = "wilcox.test", group.by = "variable", p.adjust.method = 'BH', paired = F)
res$p.signif.adj <- symnum(res$p.adj,corr = FALSE,na = FALSE,cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))


# ------------------------------------------------------------------------------------------------------------
# Extended Data Figure 9 panel b -- CD45bright subset frequency, Control vs MDS-PostAllo vs MDS-Aza (CITE-seq)
# ------------------------------------------------------------------------------------------------------------

df <- read.delim('data/EDFig9b.csv', sep = ',')
df <- df[-which(df$type == 'MDS-Aza/Ven'),]
df <- reshape2::melt(df)
df$type[which(df$type %in% c('MDS-EB1','MDS-EB2'))] <- 'MDS'
df$type <- factor(df$type, levels = c('Control','MDS','MDS-EB1','MDS-EB2','MDS-Aza','MDS-PostAllo'))
df$variable <-  gsub('\\.', ' ', df$variable)
df$variable <- factor(df$variable, levels = c('B Immature', 'B Mature', 'CD4 Naive Central Memory', 'CD4 Treg','CD4 Memory', 'CD8 Naive', 'CD8 MAIT','CD8 Memory','CD8 gdT','NK BCL11B Bright', 'NK CD56 Bright', 'NK CD56 Dim'))
df <- df[-which(df$patient %in% c('ND44_T2','MDS_2_T3','ND41_T2')),]

df$type <- factor(df$type, levels = c('Control','MDS-PostAllo','MDS','MDS-EB1','MDS-EB2','MDS-Aza'))
ggbarplot(
  df[which(df$type %in% c('Control','MDS-PostAllo', 'MDS-Aza')),], x = "variable", y = "value", #size = 0.001,
  add = c("mean_sd", 'jitter'),
  color = "type", fill = 'type',  alpha = 0.3, palette = c('#0072B2', "#fcdefc", '#a398e3'), 
  position = position_dodge(0.8)
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  ylab('Frequency') + xlab('Cell Type') + scale_color_manual(values = c('#0072B2', "#fcdefc", '#a398e3'), name = "type") + scale_y_continuous(limits = c(-0.2, 1), breaks = seq(-0.2, 1, by = 0.2))

res <- compare_means(value ~ type, data = df[which(df$type %in% c('Control','MDS-PostAllo', 'MDS-Aza')),], method = "wilcox.test", group.by = "variable", p.adjust.method = 'BH', paired = F)
res$p.signif.adj <- symnum(res$p.adj,corr = FALSE,na = FALSE,cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))


# ----------------------------------------------------------------------------------------------------
# Extended Data Figure 9 panel e -- Longitudinal TCR frequency alluvials (MDS-29, -30, -31, -33, -34)
# ----------------------------------------------------------------------------------------------------

## Single cell (CITE) alluvials
tcrs <- read.delim('data/EDFig9e.csv', sep=',')

patienttcr.rank <- function(p) {
  print(p)
  tmp <- as.data.frame(table(tcrs$tcrb[tcrs$patient == p]))
  colnames(tmp) = c('tcr', 'frequency') 
  tmp <- tmp[order(-tmp$frequency),]
  tmp$count <- tmp$frequency
  tmp$frequency <- tmp$frequency/sum(tmp$frequency)
  tmp$rank <- rank(-tmp$frequency, ties.method= "first") 
  tmp$sample <- p
  tmp$patient <- unlist(strsplit(p, '_T'))[1]
  tmp$count <- as.numeric(tmp$count)
  if(sum(tmp$count) >= 50) {
    return(tmp)
  } else {
    return()
  }
  
}

tcr_ranks <- NULL
for(samp in unique(tcrs$patient)){
  print(samp)
  tcr_ranks <- rbind(tcr_ranks, patienttcr.rank(samp))
}

library('ggalluvial')
toptcrs <- as.character(c(tcr_ranks[tcr_ranks$sample == 'ND46_T2',][1:5,]$tcr, tcr_ranks[tcr_ranks$sample == 'ND46_T4',][1:5,]$tcr, tcr_ranks[tcr_ranks$sample == 'ND46_T7',][1:5,]$tcr))
aa <- unique(toptcrs)
df <- tcr_ranks[tcr_ranks$patient == 'ND46',]
df$TCR <- as.character(df$tcr)
df$TCR[which(!(df$tcr %in% toptcrs))] <- 'Others'
df$TCR <- factor(df$TCR, levels = c(as.character(unique(toptcrs)), 'Others'))
df$sample <- factor(df$sample, levels = c('ND46_T2', 'ND46_T4', 'ND46_T7'))
ggplot(df, aes(x = sample, y = frequency, stratum = TCR, alluvium = tcr, fill = TCR)) + ggalluvial::geom_flow() + theme_classic() +
  geom_stratum(alpha = .8) + ggtitle('ND46') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(labels=c(aa, 'Others'), values = c("#cf3429", "#0072B2", "#D55E00", "#009E73", "#7a3bad", "#eb88cf", "#7bafed", "#e6b547", '#98d9cb',  "#800a01", "#043780", "#fc7905", "#035907","#480273", "#999999")) 

toptcrs <- as.character(c(tcr_ranks[tcr_ranks$sample == 'MDS_2_T1',][1:5,]$tcr, tcr_ranks[tcr_ranks$sample == 'MDS_2_T3',][1:5,]$tcr, tcr_ranks[tcr_ranks$sample == 'MDS_2_T7',][1:5,]$tcr))
aa <- unique(toptcrs)
df <- tcr_ranks[tcr_ranks$patient == 'MDS_2',]
df$TCR <- as.character(df$tcr)
df$TCR[which(!(df$tcr %in% toptcrs))] <- 'Others'
df$TCR <- factor(df$TCR, levels = c(as.character(unique(toptcrs)), 'Others'))

ggplot(df, aes(x = sample, y = frequency, stratum = TCR, alluvium = tcr, fill = TCR)) + ggalluvial::geom_flow() + theme_classic() +
  geom_stratum(alpha = .8) + ggtitle('MDS_2') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(labels=c(aa, 'Others'), values = c("#cf3429", "#0072B2", "#D55E00", "#009E73", "#7a3bad", "#eb88cf", "#800a01", "#043780", "#fc7905", "#035907","#480273","#999999"))

toptcrs <- as.character(c(tcr_ranks[tcr_ranks$sample == 'ND44_T2',][1:5,]$tcr, tcr_ranks[tcr_ranks$sample == 'ND44_T5',][1:5,]$tcr))
aa <- unique(toptcrs)
df <- tcr_ranks[tcr_ranks$patient == 'ND44',]
df$TCR <- as.character(df$tcr)
df$TCR[which(!(df$tcr %in% toptcrs))] <- 'Others'
df$TCR <- factor(df$TCR, levels = c(as.character(unique(toptcrs)), 'Others'))
df$sample <- factor(df$sample, levels = c('ND44_T2', 'ND44_T5'))
ggplot(df, aes(x = sample, y = frequency, stratum = TCR, alluvium = tcr, fill = TCR)) + ggalluvial::geom_flow() + theme_classic() +
  geom_stratum(alpha = .8) + ggtitle('ND44') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(labels=c(aa, 'Others'), values = c("#cf3429", "#0072B2", "#D55E00", "#009E73", "#7a3bad", "#eb88cf", "#7bafed", "#e6b547", '#98d9cb', "#ae92d6", "#999999")) 

toptcrs <- as.character(c(tcr_ranks[tcr_ranks$sample == 'ND41_T1',][1:5,]$tcr, tcr_ranks[tcr_ranks$sample == 'ND41_T2',][1:5,]$tcr))
aa <- unique(toptcrs)
df <- tcr_ranks[tcr_ranks$patient == 'ND41',]
df$TCR <- as.character(df$tcr)
df$TCR[which(!(df$tcr %in% toptcrs))] <- 'Others'
df$TCR <- factor(df$TCR, levels = c(as.character(unique(toptcrs)), 'Others'))
df$sample <- factor(df$sample, levels = c('ND41_T1', 'ND41_T2'))
ggplot(df, aes(x = sample, y = frequency, stratum = TCR, alluvium = tcr, fill = TCR)) + ggalluvial::geom_flow() + theme_classic() +
  geom_stratum(alpha = .8) + ggtitle('ND41') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(labels=c(aa, 'Others'), values = c("#cf3429", "#0072B2", "#D55E00", "#009E73", "#7a3bad", "#eb88cf", "#7bafed", "#e6b547", '#98d9cb', "#ae92d6", "#999999"))

toptcrs <- as.character(c(tcr_ranks[tcr_ranks$sample == 'ND92_T1',][1:5,]$tcr, tcr_ranks[tcr_ranks$sample == 'ND92_T3',][1:5,]$tcr))
aa <- unique(toptcrs)
df <- tcr_ranks[tcr_ranks$patient == 'ND92',]
df$TCR <- as.character(df$tcr)
df$TCR[which(!(df$tcr %in% toptcrs))] <- 'Others'
df$TCR <- factor(df$TCR, levels = c(as.character(unique(toptcrs)), 'Others'))
df$sample <- factor(df$sample, levels = c('ND92_T1', 'ND92_T3'))
ggplot(df, aes(x = sample, y = frequency, stratum = TCR, alluvium = tcr, fill = TCR)) + ggalluvial::geom_flow() + theme_classic() +
  geom_stratum(alpha = .8) + ggtitle('ND92') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(labels=c(aa, 'Others'), values = c("#cf3429", "#0072B2", "#D55E00", "#009E73", "#7a3bad", "#eb88cf", "#7bafed", "#e6b547", '#98d9cb', "#ae92d6", "#999999")) 



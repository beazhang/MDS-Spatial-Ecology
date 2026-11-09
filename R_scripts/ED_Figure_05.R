# ============================================================================
# Extended Data Figure 5 -- R panels
#
# Run source("_setup.R") first for the shared library() calls.
# ============================================================================


# --------------------------------------------------------------------------
# Extended Data Figure 5 panel a -- Immune aggregate size vs -PC1
# and PC1 eigenvector component calculation
# --------------------------------------------------------------------------

df <- read.delim('data/EDFig5a.csv.gz', sep = ',')
df$cat <- 'MDS'
df$cat[which(df$type %in% c('Control-Donor','Control-Lymph'))] <- 'Control'
df$cat[which(df$type %in% c('Control-TP'))] <- 'Post Allo-HCT'
df <- df[which(df$from %in% c('B Mature','B Immature','pDC','NK', 'CD8 Memory','CD8 Naive','CD4 Memory','CD4 Naive/Central Memory','CD4 Treg')),]
df1 <- dcast(data = df,formula = tls~from+neighbor,value.var = "enrichment")
df1[is.na(df1)] <- 0
library(ggfortify)
pca_res <- prcomp(df1[, 2:ncol(df1)], scale = TRUE)
df1$type <- df$cat[match(df1$tls, df$tls)]
df1$size <- df$size[match(df1$tls, df$tls)]

#PC1 eigenvector
eigenvector <- as.data.frame(-pca_res$rotation[,1])
colnames(eigenvector) <- 'rotation'
eigenvector$from <- unlist(lapply(as.character(rownames(eigenvector)), function(x) unlist(strsplit(x, '_'))[1]))
eigenvector$to <- unlist(lapply(as.character(rownames(eigenvector)), function(x) unlist(strsplit(x, '_'))[2]))
eigenvector <- dcast(eigenvector, formula = from~to, value.var = "rotation")
#write.table(eigenvector,"TLS_PC1_eigenvector_values_flipped.csv",sep=",",row.names=FALSE, quote=FALSE)


# coorelate PC1 with immune aggregate size
tmp <- data.frame('tls' = df1$tls,'pc1' = -pca_res$x[,1], 'size' = df1$size, 'type' = df1$type)
tmp$type <- factor(tmp$type, levels = c('Control','Post Allo-HCT','MDS'))
ggplot(tmp, aes(x=pc1, y=size, color = type)) + geom_point(alpha = 0.8, size = 2) +  theme_classic() + 
   scale_color_manual(values = c("#0072B2", "#c4c4c4","#BF504D"), name = "Condition") + xlab('-PC1') + ylab('Immune Aggregate Size') + 
  scale_x_continuous(limits = c(-5, 25), breaks = seq(-5, 25, by = 5)) + scale_y_continuous(limits = c(0, 1200), breaks = seq(0, 1200, by = 200))
cor.test(tmp$pc1, tmp$size)


# ------------------------------------------------------------------------------------
# Extended Data Figure 5 panel g -- MSC frequency near vs not near immune aggregate
# ------------------------------------------------------------------------------------

### Control
tls <- read.delim('data/EDFig5g_Control.csv', sep = ',')
tls <- melt(tls)

ggpaired(
  tls[which(tls$variable %in% c('Mesenchymal.Stromal')),], 
  x = "near_TLS", y = "value", alpha = 0.3,
  color = "near_TLS", palette = c("#807F7F", "#BF504D"),alpha = 0.3,
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  stat_compare_means(aes(group = near_TLS),hide.ns=T, paired=T) + ylab('Frequency within Sample') + xlab('Near Immune Aggregate') + 
  scale_color_manual(values = c("#807F7F", "#BF504D"), name = "Near Immune Aggregate") + scale_fill_manual(values = c("#807F7F", "#BF504D"), name = "Near Immune Aggregate") + scale_y_continuous(limits = c(0.0, 0.4), breaks = seq(0.0, 0.4, by = 0.1))

### MDS
tls <- read.delim('data/EDFig5g_MDS.csv', sep = ',')
tls <- melt(tls)

ggpaired(
  tls[which(tls$variable %in% c('Mesenchymal.Stromal')),], 
  x = "near_TLS", y = "value", alpha = 0.3,
  color = "near_TLS", palette = c("#807F7F", "#BF504D"),alpha = 0.3,
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  stat_compare_means(aes(group = near_TLS),hide.ns=T, paired=T) + ylab('Frequency within Sample') + xlab('Near Immune Aggregate') + 
  scale_color_manual(values = c("#807F7F", "#BF504D"), name = "Near Immune Aggregate") + scale_fill_manual(values = c("#807F7F", "#BF504D"), name = "Near Immune Aggregate") + scale_y_continuous(limits = c(0.0, 0.4), breaks = seq(0.0, 0.4, by = 0.1))

                                
# --------------------------------------------------------------------------------------
# Extended Data Figure 5 panel h -- T cell subset frequency in vs out immune aggregates
# --------------------------------------------------------------------------------------

### Control
tmp <- read.delim('data/EDFig5h_Control.csv', sep=',')
tmp <- reshape2::melt(tmp)
tmp$variable <- factor(tmp$variable, levels = c('CD8.Naive','CD8.Memory','CD4.Naive.Central.Memory','CD4.Memory','CD4.Treg'))
ggbarplot(
  tmp, x = "variable", y = "value", #size = 2,
  add = c("mean_sd", "jitter"),
  color = "in_tls", fill = 'in_tls', palette = c("#807F7F", "#BF504D"), alpha = 0.3,
  position = position_dodge(0.8)
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  ylab('Frequency within T Cells') + xlab('Cell Type') + scale_color_manual(values = c("#807F7F", "#BF504D"), name = "In Immune Aggregate") +
  scale_fill_manual(values = c("#807F7F", "#BF504D"), name = "In Immune Aggregate") + scale_y_continuous(limits = c(-0.2, 0.8), breaks = seq(-0.2, 0.8, by = 0.2))

res <- compare_means(value ~ in_tls, data = tmp, method = "wilcox.test", group.by = "variable", p.adjust.method = 'BH', paired = T)
res$p.signif.adj <- symnum(res$p.adj,corr = FALSE,na = FALSE,cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))

### MDS
tmp <- read.delim('data/EDFig5h_MDS.csv', sep=',')
tmp <- reshape2::melt(tmp)
tmp$variable <- factor(tmp$variable, levels = c('CD8.Naive','CD8.Memory','CD4.Naive.Central.Memory','CD4.Memory','CD4.Treg'))
ggbarplot(
  tmp, x = "variable", y = "value", #size = 2,
  add = c("mean_sd", "jitter"),
  color = "in_tls", fill = 'in_tls', palette = c("#807F7F", "#BF504D"), alpha = 0.3,
  position = position_dodge(0.8)
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) + 
  ylab('Frequency within T Cells') + xlab('Cell Type') + scale_color_manual(values = c("#807F7F", "#BF504D"), name = "In Immune Aggregate") +
  scale_fill_manual(values = c("#807F7F", "#BF504D"), name = "In Immune Aggregate") + scale_y_continuous(limits = c(-0.2, 0.8), breaks = seq(-0.2, 0.8, by = 0.2))

res <- compare_means(value ~ in_tls, data = tmp, method = "wilcox.test", group.by = "variable", p.adjust.method = 'BH', paired = T)
res$p.signif.adj <- symnum(res$p.adj,corr = FALSE,na = FALSE,cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns"))


# --------------------------------------------------------------------------------------------------
# Extended Data Figure 5 panel i -- Volcano: Memory T cells in MDS, in vs out immune aggregates
# --------------------------------------------------------------------------------------------------

df <- read.delim('data/EDFig5i.csv', sep = ',')
library('ggrepel')
df$diffacessible <- 'Not Significant'
df$diffacessible[df$log2FoldChange > 0.25 & df$padj < 0.05] <- paste0("Up in Aggregate")
df$diffacessible[df$log2FoldChange < -0.25 & df$padj < 0.05] <- paste0("Up outside of Aggregate")

df$DAlabel <- df$X 
df$DAlabel[df$diffacessible == 'Not Significant'] <- NA
df$DAlabel[abs(df$log2FoldChange) < 0.25] <- NA
df$DAlabel[which(!(df$DAlabel %in% c('CD244','KLRC1','KLRB1','IL7R', 'TIGIT','SELL','TCF7','CXCR4', 'CD27','CXCL9','CCL19','CXCL10','CXCL13', 'LEF1','KLRG1','PRF1','KLRD1','GZMA', 'GNLY', 'CTLA4','NKG7','CCL5','KLRD1','GZMK', 'IRF8', 'CCR7', 'CD8A','CD4', 'FOXP3')))] <- NA


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
  scale_color_manual(name = "DE Isoforms", values = mycolors1) + theme(plot.title = element_text(hjust = 0.5))  + labs(title=paste0('Memory T Cells in MDS')) + 
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 5)) + scale_x_continuous(limits = c(-1.5, 1), breaks = seq(-1.5, 1, by = 0.5))

# ============================================================================
# Extended Data Figure 1 -- R panels
#
# Run source("_setup.R") first for the shared library() calls.
# ============================================================================


# --------------------------------------------------------------------------
# Extended Data Figure 1 panel a -- Mutation oncoprint (samples x genes)
# --------------------------------------------------------------------------

library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)   # column_to_rownames; was library(tidyverse), narrowed to what is used

df <- read.delim('data/EDFig1a.csv', sep = ',')
df <- df[1:41,]

# hierarchical clustering on samples (columns in final plot)
mat <- df %>% column_to_rownames("Manuscript.ID") %>%
  mutate(across(everything(), ~ ifelse(. == "Yes", 1, 0)))
sample_order <- hclust(dist(mat, method = "binary"))$order
sample_levels <- df$Manuscript.ID[sample_order]

gene_order <- hclust(dist(t(mat), method = "binary"))$order
gene_levels <- colnames(mat)[gene_order]

# pivot to long format with clustered factor levels
df_long <- df %>%
  pivot_longer(cols = -`Manuscript.ID`, names_to = "Gene", values_to = "Mutated") %>%
  mutate(
    Manuscript.ID = factor(Manuscript.ID, levels = sample_levels),
    Gene = factor(Gene, levels = sort(gene_levels)) 
  )

# plot
ggplot(df_long, aes(x = Manuscript.ID, y = Gene, fill = Mutated)) +
  geom_tile(color = "black", lwd = 0.25) +
  scale_fill_manual(values = c("Yes" = "#405070", "No" = "#f0f0f0")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    legend.title = element_blank()
  )


# --------------------------------------------------------------------------
# Extended Data Figure 1 panel b -- CITE oncoprint
# --------------------------------------------------------------------------

df <- read.delim('data/EDFig1b.csv', sep = ',')
df <- df[,c('Patient','Type', 'Treatment', 'IPSS.M','Timepoint', 'Age','TP53.Mutant', 'SF.Mutant')]
df$tmp <- df$Timepoint
colnames(df) <-  c('Patient','Type', 'Treatment', 'IPSS-M', 'Timepoint', 'Age', 'TP53 Mutant', 'SF Mutant', 'tmp')
df <- df[,rev(c('Patient','Type', 'Treatment', 'IPSS-M', 'Timepoint','Age','TP53 Mutant', 'SF Mutant', 'tmp'))]
df[,colnames(df)] <- lapply(df[,colnames(df)] , factor)

df <- melt(df, id.var = c('Patient', 'tmp')) 

df$Patient <- factor(df$Patient, levels = c('HSS_78','HSS_79','HSS_91','HSS_179','HSS_214', 'EDTA_016','EDTA_012','EDTA_018','ND83','MDS_10',
                                            'MDS_5','MDS_9','ND150','MDS_26','ND32','ND54','MDS_12','MDS_19','ND143','ND158','MDS_18','ND166','MDS_1','MDS_30',
                                            'MDS_3','EDTA_005','MDS_27','ND120', 'EDTA_006','ND126','ND72','MDS_20',
                                            'MDS_15', 'MDS_2','ND46','ND41','MDS_28','ND92','ND44','MDS_13', 'MDS_29','MDS_16','MDS_23'))

df$value <- factor(df$value, levels = unique(df$value))
df <- df[-which(df$variable == 'Timepoint'),]
df$value <- factor(df$value, levels = c('Control','Remission','MDS','MDS-Treated','MDS-EB2','Untreated','Azacitidine','Azacitidine + Venetoclax', 'Post Allo-HCT',
                                        'Very Low','Low','Moderate Low',
                                        'Moderate High', 'High', 'Very High', 'NA','Y','N', 'Under 40','Over 60'))

ggplot(df, aes(x = tmp, y = variable, fill = value)) + 
  geom_tile(color = "black", lwd = 0.5) + 
  scale_fill_manual(values = c("N" = '#d1d1d1', "Y" = "#3b5173", 'Remission' = '#a6d4f5',
                               'MDS' = "#BD0026",'MDS-Treated' = "#e08b89",'MDS-EB2' = '#850310','Control' = '#0072B2',
                               'Under 40'= '#babcbf','40 to 60' = '#9aa0a6','Over 60' = '#404142',
                               'Untreated' = "#7258b4",'Azacitidine' = '#a398e3', 'Azacitidine + Venetoclax' = "#dbaaec", 'Post Allo-HCT' = "#fcdefc",
                               'Very Low' = "#FFFFB2", 'Low' = "#FED976", 'Moderate Low' = "#FEB24C",
                               'Moderate High' = "#FD8D3C", 'High' ="#F03B20", 'Very High' = "#BD0026")) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_grid(~Patient, drop = TRUE, scales = "free", space = "free") + 
  theme_classic() + guides(fill=guide_legend(ncol=2)) +
  theme(
    axis.title.x=element_blank(),
    panel.spacing = unit(.05, "lines"),
    panel.border = element_blank(), 
    strip.background = element_blank(),
    strip.text = element_text(angle = 90, hjust = 0)) 


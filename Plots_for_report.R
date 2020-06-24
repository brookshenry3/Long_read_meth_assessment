#Loading libraries needed

library(GenomicRanges)
library(annotatr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

#Loading shared loci csv files generated using shared_loci.R

colo <- read.csv(file = '/Users/brookshenry/Desktop/QIMR Internship/scripts/data/colo829_shared_loci.csv', sep = ',', header = TRUE)

colobl <- read.csv(file = '/Users/brookshenry/Desktop/QIMR Internship/scripts/data/colo829_bl_shared_loci.csv', sep = ',', header = TRUE)

hcc <- read.csv(file = '/Users/brookshenry/Desktop/QIMR Internship/scripts/data/hcc1937_shared_loci.csv', sep = ',', header = TRUE)

hccbl <- read.csv(file = '/Users/brookshenry/Desktop/QIMR Internship/scripts/data/hcc1937_bl_shared_loci.csv', sep = ',', header = TRUE)

pahtm <- read.csv(file = '/Users/brookshenry/Desktop/QIMR Internship/scripts/data/PAH-105TM_shared_loci.csv', sep = ',', header = TRUE)

pahbc <- read.csv(file = '/Users/brookshenry/Desktop/QIMR Internship/scripts/data/PAH-105BC_shared_loci.csv', sep = ',', header = TRUE)

#Giving all dfs a 'sample' column
colo$sample <- 'COLO829'
colobl$sample <- 'COLO829_BL'
hcc$sample <- 'HCC1937'
hccbl$sample <- 'HCC1937_BL'
pahbc$sample <- 'PAH-105BC'
pahtm$sample <- 'PAH-105TM'

#Filtering all shared loci to only include ones where the PromethION depth was > 10
colo <- subset(colo, colo$prom_called_sites > 10)
colobl <- subset(colobl, colobl$prom_called_sites > 10)
hcc <- subset(hcc, hcc$prom_called_sites > 10)
hccbl <- subset(hccbl, hccbl$prom_called_sites > 10)
pahtm <- subset(pahtm, pahtm$prom_called_sites > 10)
pahbc <- subset(pahbc, pahbc$prom_called_sites > 10)

#Resetting row names
rownames(colo) <- NULL
rownames(colobl) <- NULL
rownames(hcc) <- NULL
rownames(hccbl) <- NULL
rownames(pahtm) <- NULL
rownames(pahbc) <- NULL


###################

#Making one large shared loci df
shared.loci <- rbind(colo, colobl, hcc, hccbl, pahbc, pahtm)

#Filtering PromthION positions where there are fewer than 10 reads present - redundant step
shared.loci <- subset(shared.loci, shared.loci$prom_called_sites > 10)

#rm(colobl, hcc, hccbl, pahbc, pahtm)

#Cleaning things up, factorizing sample
shared.loci <- shared.loci[, c(1, 2, 5, 9, 10)]
shared.loci$sample <- as.factor(shared.loci$sample)

##################

#Correlation testing, 
colo.corr <- cor.test(colo$epic_beta, colo$prom_meth_freq)
colobl.corr <- cor.test(colobl$epic_beta, colobl$prom_meth_freq)
hcc.corr <- cor.test(hcc$epic_beta, hcc$prom_meth_freq)
hccbl.corr <- cor.test(hccbl$epic_beta, hccbl$prom_meth_freq)
pahbc.corr <- cor.test(pahbc$epic_beta, pahbc$prom_meth_freq)
pahtm.corr <- cor.test(pahtm$epic_beta, pahtm$prom_meth_freq)

#Plotting correlation between EPIC and PromethION methylation values

#c.test <- cor.test(shared.loci$epic_beta, shared.loci$prom_meth_freq)

#c <- cor(shared.loci$epic_beta, shared.loci$prom_meth_freq)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

#title <- sprintf("N = %d r = 0.9500", nrow(shared.loci), c)

ggplot(shared.loci, aes(epic_beta, prom_meth_freq)) +
  geom_bin2d(bins=25) + scale_fill_gradientn(colors=r, trans="log10", name = 'Count') +
  xlab("EPIC Methylation Frequency") +
  ylab("PromethION Methylation Frequency") +
  theme_bw(base_size=20) +
  #ggtitle(title)
  facet_wrap(~sample) +
  theme(axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))

rm(rf, r, c.test)

##################

#Plotting distributions of EPIC and PromethION methylation values 

dist <- data.frame(value = c(colo$epic_beta, colo$prom_meth_freq, colobl$epic_beta, colobl$prom_meth_freq, hcc$epic_beta, hcc$prom_meth_freq, hccbl$epic_beta, hccbl$prom_meth_freq, pahbc$epic_beta, pahbc$prom_meth_freq, pahtm$epic_beta, pahtm$prom_meth_freq),
                   method = c(rep('EPIC', nrow(colo)), rep('PromethION', nrow(colo)), rep('EPIC', nrow(colobl)), rep('PromethION', nrow(colobl)), rep('EPIC', nrow(hcc)), rep('PromethION', nrow(hcc)), rep('EPIC', nrow(hccbl)), rep('PromethION', nrow(hccbl)), rep('EPIC', nrow(pahbc)), rep('PromethION', nrow(pahbc)), rep('EPIC', nrow(pahtm)), rep('PromethION', nrow(pahtm))),
                   sample = c(colo$sample, colobl$sample, hcc$sample, hccbl$sample, pahbc$sample, pahtm$sample))

ggplot(data=dist, aes(x=value, group=method, fill=method)) +
  geom_density(adjust=1.5, alpha=.4) + 
  #ggtitle('EPIC vs. PromethION values') +
  xlab('Methylation Value') +
  ylab('Density') +
  guides(fill=guide_legend(title=NULL)) +
  facet_wrap(~sample)
 
#Testing difference in methylation value distributions with Kolmogorov-Smirnov (KS) tests
ks.test(colo$epic_beta, colo$prom_meth_freq)
ks.test(colobl$epic_beta, colobl$prom_meth_freq)
ks.test(hcc$epic_beta, hcc$prom_meth_freq)
ks.test(hccbl$epic_beta, hccbl$prom_meth_freq)
ks.test(pahbc$epic_beta, pahbc$prom_meth_freq)
ks.test(pahtm$epic_beta, pahtm$prom_meth_freq)


rm(dist)

##################

#Annotating shared loci with genomic annotations for annotation enrichment analysis

annots <- c('hg38_basicgenes', 'hg38_genes_intergenic', 'hg38_cpg_islands', 'hg38_cpg_inter', 'hg38_cpgs', 'hg38_enhancers_fantom')

annotations = build_annotations(genome = 'hg38', annotations = annots)


#Converting shared loci to GRanges for annotation
colo.gr <- makeGRangesFromDataFrame(colo,
                                      keep.extra.columns = TRUE,
                                      ignore.strand = TRUE,
                                      seqnames.field = 'chromosome',
                                      start.field = 'position',
                                      end.field = 'position')

colobl.gr <- makeGRangesFromDataFrame(colobl,
                                           keep.extra.columns = TRUE,
                                           ignore.strand = TRUE,
                                           seqnames.field = 'chromosome',
                                           start.field = 'position',
                                           end.field = 'position')

hcc.gr <- makeGRangesFromDataFrame(hcc,
                                      keep.extra.columns = TRUE,
                                      ignore.strand = TRUE,
                                      seqnames.field = 'chromosome',
                                      start.field = 'position',
                                      end.field = 'position')

hccbl.gr <- makeGRangesFromDataFrame(hccbl,
                                      keep.extra.columns = TRUE,
                                      ignore.strand = TRUE,
                                      seqnames.field = 'chromosome',
                                      start.field = 'position',
                                      end.field = 'position')

pahbc.gr <- makeGRangesFromDataFrame(pahbc,
                                      keep.extra.columns = TRUE,
                                      ignore.strand = TRUE,
                                      seqnames.field = 'chromosome',
                                      start.field = 'position',
                                      end.field = 'position')

pahtm.gr <- makeGRangesFromDataFrame(pahtm,
                                      keep.extra.columns = TRUE,
                                      ignore.strand = TRUE,
                                      seqnames.field = 'chromosome',
                                      start.field = 'position',
                                      end.field = 'position')

#Annotating shared loci GRanges objects with genomic annotations
colo.annot <- annotate_regions(regions = colo.gr,
                                 annotations = annotations,
                                 ignore.strand = TRUE,
                                 quiet = FALSE)

colobl.annot <- annotate_regions(regions = colobl.gr,
                             annotations = annotations,
                             ignore.strand = TRUE,
                             quiet = FALSE)

hcc.annot <- annotate_regions(regions = hcc.gr,
                                 annotations = annotations,
                                 ignore.strand = TRUE,
                                 quiet = FALSE)

hccbl.annot <- annotate_regions(regions = hccbl.gr,
                                 annotations = annotations,
                                 ignore.strand = TRUE,
                                 quiet = FALSE)

pahbc.annot <- annotate_regions(regions = pahbc.gr,
                                 annotations = annotations,
                                 ignore.strand = TRUE,
                                 quiet = FALSE)

pahtm.annot <- annotate_regions(regions = pahtm.gr,
                                 annotations = annotations,
                                 ignore.strand = TRUE,
                                 quiet = FALSE)

#Converting annotated GRanges objects back to data frames
colo.annot.df <- as.data.frame(colo.annot)
colobl.annot.df <- as.data.frame(colobl.annot)
hcc.annot.df <- as.data.frame(hcc.annot)
hccbl.annot.df <- as.data.frame(hccbl.annot)
pahbc.annot.df <- as.data.frame(pahbc.annot)
pahtm.annot.df <- as.data.frame(pahtm.annot)

#Making annotations factors
colo.annot.df$annot.type <- as.factor(colo.annot.df$annot.type)
colobl.annot.df$annot.type <- as.factor(colobl.annot.df$annot.type)
hcc.annot.df$annot.type <- as.factor(hcc.annot.df$annot.type)
hccbl.annot.df$annot.type <- as.factor(hccbl.annot.df$annot.type)
pahbc.annot.df$annot.type <- as.factor(pahbc.annot.df$annot.type)
pahtm.annot.df$annot.type <- as.factor(pahtm.annot.df$annot.type)

#Getting difference in methylation calls between EPIC and PromethION values
colo.annot.df$meth_diff <- colo.annot.df$epic_beta - colo.annot.df$prom_meth_freq
colobl.annot.df$meth_diff <- colobl.annot.df$epic_beta - colobl.annot.df$prom_meth_freq
hcc.annot.df$meth_diff <- hcc.annot.df$epic_beta - hcc.annot.df$prom_meth_freq
hccbl.annot.df$meth_diff <- hccbl.annot.df$epic_beta - hccbl.annot.df$prom_meth_freq
pahbc.annot.df$meth_diff <- pahbc.annot.df$epic_beta - pahbc.annot.df$prom_meth_freq
pahtm.annot.df$meth_diff <- pahtm.annot.df$epic_beta - pahtm.annot.df$prom_meth_freq

d <- 0.25 #Setting threshold to consider methylation calls as non-concordant

#Finding shared loci that fall below the threshold set above
colo.diff.meth <- subset(colo.annot.df, colo.annot.df$meth_diff > d | colo.annot.df$meth_diff < -(d))
colobl.diff.meth <- subset(colobl.annot.df, colobl.annot.df$meth_diff > d | colobl.annot.df$meth_diff < -(d))
hcc.diff.meth <- subset(hcc.annot.df, hcc.annot.df$meth_diff > d | hcc.annot.df$meth_diff < -(d))
hccbl.diff.meth <- subset(hccbl.annot.df, hccbl.annot.df$meth_diff > d | hccbl.annot.df$meth_diff < -(d))
pahbc.diff.meth <- subset(pahbc.annot.df, pahbc.annot.df$meth_diff > d | pahbc.annot.df$meth_diff < -(d))
pahtm.diff.meth <- subset(pahtm.annot.df, pahtm.annot.df$meth_diff > d | pahtm.annot.df$meth_diff < -(d))

#Summarizing all annotation types in both all shared loci and those identified above as having non-concordant methylation calls
colo.annot.sum <- summary(colo.annot.df$annot.type)
colo.dm.sum <- summary(colo.diff.meth$annot.type)

colobl.annot.sum <- summary(colobl.annot.df$annot.type)
colobl.dm.sum <- summary(colobl.diff.meth$annot.type)

hcc.annot.sum <- summary(hcc.annot.df$annot.type)
hcc.dm.sum <- summary(hcc.diff.meth$annot.type)

hccbl.annot.sum <- summary(hccbl.annot.df$annot.type)
hccbl.dm.sum <- summary(hccbl.diff.meth$annot.type)

pahbc.annot.sum <- summary(pahbc.annot.df$annot.type)
pahbc.dm.sum <- summary(pahbc.diff.meth$annot.type)

pahtm.annot.sum <- summary(pahtm.annot.df$annot.type)
pahtm.dm.sum <- summary(pahtm.diff.meth$annot.type)

#Getting summary of annotation types ready for bar plots
diff.annot.enrich <- data.frame("annot" = rep(c('CpG_inter', 'CpG_islands', 'CpG_shelves', 'CpG_shores', 'FANTOM_enhnacers', '1to5kb', '3UTR', '5UTR', 'exons', 'intergenic', 'introns', 'promoters'), 6),
                                "group" = c(rep("total", 12), rep("diff_meth", 12), rep("total", 12), rep("diff_meth", 12), rep("total", 12), rep("diff_meth", 12), rep("total", 12), rep("diff_meth", 12), rep("total", 12), rep("diff_meth", 12), rep("total", 12), rep("diff_meth", 12)),
                                "number_annots" = c(colo.annot.sum, colo.dm.sum, colobl.annot.sum, colobl.dm.sum, hcc.annot.sum, hcc.dm.sum, hccbl.annot.sum, hccbl.dm.sum, pahbc.annot.sum, pahbc.dm.sum, pahtm.annot.sum, pahtm.dm.sum),
                                "total" = c(rep(nrow(colo.annot.df), 12), rep(nrow(colo.diff.meth), 12), rep(nrow(colobl.annot.df), 12), rep(nrow(colobl.diff.meth), 12), rep(nrow(hcc.annot.df), 12), rep(nrow(hcc.diff.meth), 12), rep(nrow(hccbl.annot.df), 12), rep(nrow(hccbl.diff.meth), 12), rep(nrow(pahbc.annot.df), 12), rep(nrow(pahbc.diff.meth), 12), rep(nrow(pahtm.annot.df), 12), rep(nrow(pahtm.diff.meth), 12)),
                                "sample" = c(rep('COLO829', 24), rep('COLO829_BL', 24), rep('HCC1937', 24), rep('HCC1937_BL', 24), rep('PAH-105BC', 24), rep('PAH-105TM', 24)),
                                "percent" = NA)

diff.annot.enrich$percent <- diff.annot.enrich$number_annots / diff.annot.enrich$total * 100

#Plotting annotation enrichment
ggplot(diff.annot.enrich, aes(fill=group, y=percent, x=annot)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x='Annotation', y='Percent', fill='Group') + 
  theme(legend.position = 'none') +
  #ggtitle("Annotation enrichment in loci with a methylation difference of |>0.75| between EPIC and PromethION") +
  facet_wrap(~sample)


#####

#Plotting distribution of differences in methylation values - used to set threshold above for non-concordant methylation calls

colo$meth_diff <- colo$epic_beta - colo$prom_meth_freq
colobl$meth_diff <- colobl$epic_beta - colobl$prom_meth_freq
hcc$meth_diff <- hcc$epic_beta - hcc$prom_meth_freq
hccbl$meth_diff <- hccbl$epic_beta - hccbl$prom_meth_freq
pahbc$meth_diff <- pahbc$epic_beta - pahbc$prom_meth_freq
pahtm$meth_diff <- pahtm$epic_beta - pahtm$prom_meth_freq


dist <- data.frame(value = c(colo$meth_diff, colobl$meth_diff, hcc$meth_diff, hccbl$meth_diff, pahbc$meth_diff, pahtm$meth_diff),
                   sample = c(colo$sample, colobl$sample, hcc$sample, hccbl$sample, pahbc$sample, pahtm$sample))

ggplot(data=dist, aes(x=value, fill=sample)) +
  geom_density(adjust=1.5, alpha=.4) + 
  #ggtitle('EPIC vs. PromethION values') +
  xlab('Difference in Methylation Call') +
  ylab('Density') +
  guides(fill=guide_legend(title=NULL)) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype='dashed')

rm(dist)



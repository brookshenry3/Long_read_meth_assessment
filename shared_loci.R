#This is the full workflow for the methylation data from both the EPIC array and PromethION 

####Packages####
library(GenomicRanges)
library(annotatr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

#####Loading/formatting data####

load(file = 'path/to/epic_meth_w_betas_gr') #Loading the EPIC GRanges object generated in minfi.R
epic.meth <- as.data.frame(epic.meth.gr) #Converting the granges object above into a df

p.meth <- read.table(file = 'path/to/longread/meth_frequencies.tsv', sep = '\t', header = TRUE) #This tsv file loaded here was generated using f5c's meth freq function

#Cleaning/reformatting data


p.meth$start <- p.meth$end #ok for some reason this tsv has no "start' positions - they're listed as 0. So I'll duplicate the end positions to get the proper formattings. 

#Converting the PromethION data to a GRanges object
p.meth.gr <- makeGRangesFromDataFrame(p.meth,
                                      keep.extra.columns = TRUE,
                                      ignore.strand = TRUE,
                                      seqnames.field = 'chromosome',
                                      start.field = 'start',
                                      end.field = 'end')

epic.meth.gr <- sortSeqlevels(epic.meth.gr) #Sorting the chromosomes so that they're in the proper order
p.meth.gr <- sortSeqlevels(p.meth.gr) #Sorting the chromosomes so that they're in the proper order

epic.meth.gr <- sort(epic.meth.gr)
p.meth.gr <- sort(p.meth.gr)

shifted.p <- GenomicRanges::shift(p.meth.gr, 1) #Shifting the PromethION data by 1bp so that it lines up with the EPIC data

####Finding common loci between the two methods####

#Now to find overlapping ranges

common.loci.e <- subsetByOverlaps(epic.meth.gr, shifted.p, type = 'equal')
common.loci.p <- subsetByOverlaps(shifted.p, common.loci.e, type = 'equal')

cl.e.df <- as.data.frame(common.loci.e)
cl.p.df <- as.data.frame(common.loci.p)

#So weirdly above even after doing subsetting there is a 14 row difference between the two of them

pos.e <- cl.e.df[, c(1:3)]
pos.p <- cl.p.df[, c(1:3)]

#diff <- duplicated(pos.e)

pos.e <- pos.e[!duplicated(pos.e), ]
#rownames(pos.e) <- c(1:nrow(pos.e))

cl.e.df <- subset(cl.e.df, rownames(cl.e.df) %in% rownames(pos.e))

#Creating/saving a shared loci dataframe 
shared.loci <- data.frame(chromosome = cl.e.df$seqnames,
                          position = cl.e.df$start,
                          epic_probe = cl.e.df$Name,
                          epic_seq = cl.e.df$ProbeSeqA,
                          epic_beta = cl.e.df$X204228990124_R06C01,
                          prom_seq = cl.p.df$group_sequence,
                          prom_called_sites = cl.p.df$called_sites,
                          prom_meth_sites = cl.p.df$called_sites_methylated,
                          prom_meth_freq = cl.p.df$methylated_frequency)

write.csv(shared.loci, "path/to/write/shared_loci.csv", row.names = FALSE)








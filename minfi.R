#Packages used by minfi
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(GenomicRanges)

#Packages used to facilitate liftover
library(rtracklayer)

#Loading in the raw green/red .idat files for the EPIC array
rgset <- read.metharray(basenames = 'Path/to/idat/files')

qcReport(rgset, pdf= 'path/to/output/qcreport.pdf')

####Preprocessing the data (Illumina method)####

mset.illumina <- preprocessIllumina(rgset, bg.correct = TRUE, normalize = 'controls')

#Now with that preprocessing done I can convert the mset.illumina into a GenomicRatioSet (also getting the betas here)

rset.illumina <- ratioConvert(mset.illumina)

#betas.illumina <- getBeta(rset.illumina)

grset.illumina <- mapToGenome(rset.illumina, mergeManifest = TRUE)

#Important to note that minfi appears to be working with hg19, will need to perform a liftover to get it to hg38

#Now removing SNP probes

grset.illumina <- dropLociWithSnps(grset.illumina, snps=c('SBE','CpG'), maf=0)

betas.illumina <- getBeta(grset.illumina)

#Now converting the grset into a GRanges object

data.gr <- granges(grset.illumina)

#Now to liftover the data from hg19 to hg38

chain <- import.chain('path/to/hg19ToHg38.over.chain')

datahg38.gr <- liftOver(data.gr, chain)

#Converting the data to a dataframe to clean it up

df.epic.meth <- as.data.frame(datahg38.gr)
df.epic.meth <- df.epic.meth[ ,c(3, 4, 5, 9, 12)]
df.epic.meth$beta <- subset(betas.illumina, rownames(betas.illumina) %in% df.epic.meth$Name)

#Finally after I've added the beta values to the probe location and sequence info above I can convert the df back into a GRanges object

epic.meth.gr <- makeGRangesFromDataFrame(df.epic.meth,
                                         keep.extra.columns = TRUE,
                                         seqnames.field = 'seqnames',
                                         start.field = 'start',
                                         end.field = 'end')



#Saving all of the data
save(betas.illumina, file = 'path/to/save/data/illumina_beta_matrix')
save(data.gr, file = 'path/to/save/data/EPIC_meth_hg19')
save(datahg38.gr, file = 'path/to/save/data/EPIC_meth_hg38')
save(epic.meth.gr, file = 'path/to/save/data/epic_meth_w_betas_gr')

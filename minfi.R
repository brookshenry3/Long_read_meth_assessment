#Packages used by minfi
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(GenomicRanges)

#Packages used to facilitate liftover
library(rtracklayer)

#Loading in the raw green/red .idat files for the EPIC array from the colo829 data
rgset <- read.metharray(basenames = '/mnt/lustre/working/genomeinfo/data/20200303_AMGP/204228990124_R06C01_Grn.idat')

# #Getting only the methylated/unmethylated signals
# mset <- preprocessRaw(rgset)
# 
# #Now getting beta and m values
# rset <- ratioConvert(mset, what = 'both', keepCN = TRUE)
# 
# #Now just getting a matrix containing the beta values
# betas <- getBeta(rset)
# 
# #Running QC on the mset
# qc <- getQC(mset)
# plotQC(qc)
# densityPlot(mset)

qcReport(rgset, pdf= '/working/lab_nicw/kevinH/projects/methylation_assessment/data/minfi_meth_calls/CO9a_qcreport.pdf')

#Based on the above the data looks pretty good, not sure what QC I need to do

####Preprocessing the data (Illumina method)####

#Might as well try doing the illumina preprocessing

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

chain <- import.chain('/mnt/lustre/working/lab_nicw/kevinH/projects/methylation_assessment/general_scripts/hg19ToHg38.over.chain')

datahg38.gr <- liftOver(data.gr, chain)

#Workflow now includes the EPIC part of 'formatting_data.R' so that script is now redundant


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
#save(betas, file = '/mnt/lustre/working/lab_nicw/kevinH/projects/methylation_assessment/data/minfi_meth_calls/colo829_bl_raw_beta_matrix')
save(betas.illumina, file = '/mnt/lustre/working/lab_nicw/kevinH/projects/methylation_assessment/data/minfi_meth_calls/CO9a_illumina_beta_matrix')
save(data.gr, file = '/mnt/lustre/working/lab_nicw/kevinH/projects/methylation_assessment/data/minfi_meth_calls/CO9a_EPIC_meth_hg19')
save(datahg38.gr, file = '/mnt/lustre/working/lab_nicw/kevinH/projects/methylation_assessment/data/minfi_meth_calls/CO9a_EPIC_meth_hg38')
save(epic.meth.gr, file = '/mnt/lustre/working/lab_nicw/kevinH/projects/methylation_assessment/data/minfi_meth_calls/CO9a_epic_meth_w_betas_gr')



#load(file = '/mnt/lustre/working/lab_nicw/kevinH/projects/methylation_assessment/colo_829/colo829_EPIC_data/EPIC_meth_hg19')

#load(file = '/mnt/lustre/working/lab_nicw/kevinH/projects/methylation_assessment/colo_829/colo829_EPIC_data/EPIC_meth_hg38')


#hg19 <- as.data.frame(data.gr)
#hg38 <- as.data.frame(datahg38.gr)










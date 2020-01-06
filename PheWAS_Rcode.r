
#RUN A PHENOME-WIDE ASSOCIATION STUDY USING THE R PACKAGE, PHEWAS, FOR ARIVALE DATA (environment is my Jupyter notebook):

install.packages("devtools")
install.packages("meta")
install.packages("ggrepel")
install.packages("DT")

library(devtools)
library(meta)
library(ggrepel)
library(DT)

install.packages("PheWAS")
library(PheWAS)

#upload data file that includes client ID, clinical chemistries, proteins, metabolites, demographics, and genetic ancestry principal components), merged from Arivale Snapshots previously. Baseline blood-draw data only:
alldata <- read.csv(file="/notebooks/APOEanalyses/May2019_alldata_baseline.csv", sep=",", header=T)


#Upload your genotype data (obtained from WGS using VIPR for Arivale participants, see John Earles):
rs429358 <- read.csv(file="/notebooks/APOEanalyses/rs429358.csv", sep=",", header=T)


#create genotype file, which has two columns: client id and numeric genotype, with 0 as major allele homozygote, 1 as het, 2 as minor allele homozygote. If missingness in the genotype file, remove it first.
rs429358 <- subset(rs429358, select=c(public_client_id, genotype))
rs429358$geno <- (ifelse(rs429358$genotype=="T/T", 0,
                        ifelse(rs429358$genotype=="C/T", 1, 2)))
rs429358$genotype<-NULL


#Create covariate file with sex, age, and the first four principal components. I need two covariate files because the clinical chemistries were obtained from two different vendors, which need to be adjusted for. But the proteins and metabolites do not need to be adjusted for vendor. 

covariates1 <- subset(alldata, select=c(public_client_id, female, age, vendor_id, PC1, PC2, PC3, PC4))
covariates1B <- subset(alldata, select=c(public_client_id, female, age, PC1, PC2, PC3, PC4))

#write out covariates for use in SNP-sex interaction analyses:
write.table(covariates1, file="covers.csv", sep=",")


#log transform the phenotypes to correct for right skew and minimize outliers. 
#first separate the proteomics from the chems & mets, proteins were already transformed and batch-corrected by Olink:

public_client_id <- alldata2[,1]
head(public_client_id)

chems <- subset(alldata1, select=c(24:85))
logchems <- log(chems+1)	#add a 1 to deal with zeroes
logchems <- cbind(public_client_id, logchems)
head(logchems)
dim(logchems)

prots <- subset(alldata1, select=c(public_client_id, 841:2033))
head(prots)
dim(prots)

mets <- subset(alldata, select=c(88:840))
logmets <- log(mets+1)
logmets <- cbind(public_client_id, logmets)
head(logmets)
dim(logmets)

#going to run phewas separately for chemistries because of vendor effect. Merge log values of proteins & metabolites for non-vendor-adjusted phewas run
logprotmets <- merge(prots, logmets, all=TRUE)

#output the log-transformed phenotype files for use in SNP-sex interaction analysis:
write.table(logchems, file="log_chemistries.csv", sep=",")
write.table(logprotmets, file="log_protmets.csv", sep=",")


##RUN PHEWAS##

#pheinfo file made in excel, upload manually
pheinfo <- read.csv(file="pheinfo_June_2019.csv", sep=",", header=T, na.strings="", stringsAsFactors=FALSE)

#First run for chemistries, will return beta coefficient for each analyte
results<-phewas(phenotypes=logchems, genotypes=rs429358, covariates=covariates1, cores=1, significance.threshold=c("fdr"))
results_d<-addPhecodeInfo(results)
results_d[order(results_d$p),]

#Next run for proteins and metabolites
resultsB<-phewas(phenotypes=logprotmets, genotypes=s429358, covariates=covariates1B, cores=1, significance.threshold=c("fdr"))
resultsB_d<-addPhecodeInfo(resultsB)
resultsB_d[order(resultsB_d$p),]

#bind the results together so I can rerun the FDR adjustment manually with all analyses and output q-values
results_all <- rbind(results, resultsB)

adjustp <- p.adjust(results_all$p, "fdr")
adjustresults <- cbind(adjustp, results_all)

#write out results for SNP:
write.table(adjustresults, file="adjusted_rs429358_phewas.csv", sep=",")

#illustrate results with a Manhattan plot, output pdf
pdf(file="rs429358_phewas_plot.pdf", width=7, height=7)
phewasManhattan(adjustresults, annotate.level=FALSE, title="rs429358 PheWAS results", y.axis.interval=2, point.size=1, size.x.labels=10, size.y.labels=14)
dev.off()



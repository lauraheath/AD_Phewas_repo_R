#upload data files that were generated in the SNP_sex_interactions_r.txt step.
#This file has everything (log transformed analyses, PCs, covariates) together.

alldata_ln <- read.csv(file="/notebooks/baselinedata_logall.csv", sep=",", header=T)



#create sex and age labels
alldata_ln$gender <- ifelse(alldata$female==0, "men", "women")

alldata_ln$age_labels <- ifelse(alldata_ln$age_cats==0, "18-29",
                            ifelse(alldata_ln$age_cats==1, "30-39",
                                  ifelse(alldata_ln$age_cats==2, "40-49",
                                        ifelse(alldata_ln$age_cats==3, "50-59",
                                              ifelse(alldata_ln$age_cats==4, "60-69", "70+")))))





install.packages("ggplot2")
library(ggplot2)



#merge alldata file with genotype data:
alldata_ln2 <- merge(alldata, rs429358, by='public_client_id')
head(alldata_ln2)
dim(alldata_ln2)

#check the genotypes:
table(alldata_ln2$genotype)

#I want to reverse the order of the boxes so the minor allele is last instead of first
levels(alldata_ln2$genotype)
alldata_ln2$genotype <- factor(alldata_ln2$genotype, levels=rev(levels(alldata_ln2$genotype)))
levels(alldata_ln2$genotype)

#Create a graph that contains genotye-coded box plots by age group for an analyte of interest and write out a pdf of the plot:
pdf(file="boxplot_rs429358_LDLcholcalc.pdf", width=8, height=7)
p <- ggplot(alldata_ln2, aes(x=age_labels, LDL.CHOL.CALCULATION, fill=genotype)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=20)) + xlab("Age group (decade)") + ylab("LDL cholesterol")  + theme_classic(base_size=20) + scale_fill_manual(breaks = c("T/T", "C/T", "C/C"), values=c("white", "gray", "black"), name="rs429358 genotype")
p <- p + theme(axis.text.x=element_text(angle=90)) #+ ylim(9, 14)
p
dev.off()


#Two-faceted plot by gender:
pdf(file="boxplot_rs429358_LDLL_by_sex.pdf")
p <- ggplot(alldata_ln2, aes(x=genotype, LDL.CHOL.CALCULATION, fill=genotype)) 
p <- p + geom_boxplot() + theme(axis.text=element_text(size=14)) + theme_bw(base_size=18) + guides(fill=guide_legend(title="rs429358"))
p + scale_fill_manual(breaks = c("T/T", "C/T", "C/C"), values=c("white", "gray", "black")) + facet_wrap(~gender)
p
dev.off()

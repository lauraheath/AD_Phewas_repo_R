
#upload data files that were generated in the PheWAS_Rcode step. Contains covariates, log-transformed chemistries, and log-transformed proteins & metabolites files:
covars <- read.csv(file="/notebooks/covars.csv", sep=",", header=T)

logchems1 <- read.csv(file="/notebooks/APOEanalyses/May2019_logchems.csv", sep=",", header=T)
logprotmets1 <- read.csv(file="/notebooks/APOEanalyses/May2019_logprotmets.csv", sep=",", header=T)

#upload SNP data:
rs429358 <- read.csv(file="/notebooks/APOEanalyses/rs429358.csv", sep=",", header=T)
rs429358 <- subset(rs429358, select=c(public_client_id, genotype))
rs429358$geno <- (ifelse(rs429358$genotype=="T/T", 0,
                        ifelse(rs429358$genotype=="C/T", 1, 2)))



#need different adjustments for chemistries due to vendor issues, so make 2 files
logchems <- merge(logchems1, covars, all=TRUE)
logprotmets <- merge(logprotmets1, covars, all=TRUE)

#merge these files for creating box plots later, write the file out:
alldata_ln <- merge(logchems1, logprotmets1, all=TRUE)
alldata_ln <- merge(alldata_ln, covars, all=TRUE)
write.table(alldata_ln, file="baselinedata_logall.csv", sep=",")



#several chemistries were only collected with one of the vendors. take these out of logchems file and merge them with 
#proteins & metabolites file, so that the row numbers are stable due to vendor_id getting dropped in these analytes
logchems_sub <- subset(logchems, select=c(public_client_id, DPA, GFR..MDRD..AFRICAN.AM, HDL.PARTICLE.NUMBER, LDL_SIZE, LINOLEIC_ACID, OMEGA_3_TOTAL, OMEGA_6_TOTAL))
logprotmets2<-merge(logchems_sub, logprotmets)


logchems2<-logchems
logchems2$DPA<-NULL
logchems2$GFR..MDRD..AFRICAN.AM<-NULL
logchems2$HDL.PARTICLE.NUMBER<-NULL
logchems2$LDL_SIZE<-NULL
logchems2$LINOLEIC_ACID<-NULL
logchems2$OMEGA_3_TOTAL<-NULL
logchems2$OMEGA_6_TOTAL<-NULL

#output the column names in order:
chemcolnames <- colnames(logchems2)
write.table(chemcolnames, file="colnames_chems.csv", sep=",")
protmets_names <- colnames(logprotmets2)
write.table(protmets_names, file="colnames_protmets.csv", sep=",")


##RUN LINEAR REGRESSION ON ALL ANALYTES WITH AN INTERACTION TERM, CHEMISTRIES SEPARATE FROM PROTEINS/METS THEN REJOIN##

#first join the genotype data to the data files
logchems2 <- merge(logchems2, rs429358, by='public_client_id')
logprotmets2 <- merge(logprotmets2, rs429358, by='public_client_id')


#columns 2 through 56 in logchems2 has all clinical chemistries, columns 2-1954 have all proteins & metabolites in logprotmets2 file:
fits1 <- lapply(logchems2[,2:56], function(x) lm(x~rs429358*female + age + vendor_id + PC1 + PC2 + PC3 + PC4, data=logchems2))
fits2 <- lapply(logprotmets2[,2:1954], function(x) lm(x~rs429358*female + age + PC1 + PC2 + PC3 + PC4, data=logprotmets2))   

summaries1 <- lapply(fits1, summary)
summaries1 <- lapply(summaries1, function(x) x$coefficients[, c(1,4)])
    
summaries2 <- lapply(fits2, summary)
summaries2 <- lapply(summaries2, function(x) x$coefficients[, c(1,4)])
    
lm1 <- ldply(summaries1, data.frame)
#need to label the rows: first number rows within each analyte, then add labels
colnames(lm1)[1] <- "analyte"

lm2 <- ldply(summaries2, data.frame)
colnames(lm2)[1] <- "analyte"

lm1$rownum <- sequence(rle(lm1$analyte)$lengths)
lm2$rownum <- sequence(rle(lm2$analyte)$lengths)

#here are the labels
lm1$coef <- ifelse(lm1$rownum==1, "intercept",
                  ifelse(lm1$rownum==2, "geno",
                               ifelse(lm1$rownum==3, "female",
                                      ifelse(lm1$rownum==4, "age", 
                                             ifelse(lm1$rownum==5, "vendor_id", 
                                                    ifelse(lm1$rownum==6, "PC1", 
                                                           ifelse(lm1$rownum==7, "PC2", 
                                                                  ifelse(lm1$rownum==8, "PC3",
                                                                         ifelse(lm1$rownum==9, "PC4", "geno:female")))))))))
lm2$coef <- ifelse(lm2$rownum==1, "intercept",
                  ifelse(lm2$rownum==2, "geno",
                               ifelse(lm2$rownum==3, "female",
                                      ifelse(lm2$rownum==4, "age", 
                                             ifelse(lm2$rownum==5, "PC1", 
                                                    ifelse(lm2$rownum==6, "PC2", 
                                                           ifelse(lm2$rownum==7, "PC3",
                                                                  ifelse(lm2$rownum==8, "PC4", "geno:female"))))))))
                                                                                          
                                                                                                  


#calculate fdr adjusted p-values for the interaction coefficient only:
#first bind the chems & prot/met dataframes together:
lm_all <- rbind(lm1, lm2)
dim(lm_all)
fdr_int <- subset(lm_all, lm_all$coef=="geno:female")
head(fdr_int)
dim(fdr_int)

adjustp_int <- p.adjust(fdr_int$Pr...t.., "fdr")
fdrP <- cbind(adjustp_int, fdr_int)

fdrP$Estimate<-NULL
fdrP$Pr...t..<-NULL
fdrP$rownum<-NULL

lm_adjusted <- merge(lm_all, fdrP, all=TRUE)
#write out all results:
write.table(lm_adjusted, file="rs429358_sex_interaction.csv", sep=",")


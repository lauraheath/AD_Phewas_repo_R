PheWAS_Rcode.r
Perform phenome-wide association study using R package "PheWAS" with Arivale data Running the R package “PheWAS”
Directions from Vanderbilt: https://www.vumc.org/cpm/cpm-blog/phewas-r-package
Github: https://github.com/PheWAS/PheWAS
Paper to be cited: https://www.nature.com/articles/nbt.2749
This package performs a phenome-wide association study, automatically performs linear or logistic regression based on input, returns main effect beta coefficients/ORs, SEs, p-values, allele frequencies, HWE, corrects for multiple tests via multiple methods, and creates attractive Manhattan plots.
Input requires:
	1.	Genotype file, containing subject ID and genotype, coded numerically 0,1,2 for typical SNP analyses (but this value can be any predictor).
	2.	Covariate file, containing subject ID and any covariates (sex, age, BMI, comoborbidities, medication use, etc).
	3.	Phenotype file, containing subject ID and each outcome (ICD-9 code status or quantitative analyte).
	4.	PheCode file: there is a file that comes with the PheWAS package to translate ICD9-codes into outcomes for the analysis. This file also makes it possible to create pretty Manhattan plots of results. If working with outcomes other than ICD9 codes, then make and upload your own file, labeled “pheinfo,” and it will be incorporated into the PheWAS analysis. For Arivale data, I created this file in excel. Requires five columns: Phecode = ICD9 code or analyte name, exactly as it is in the phenotype file. Description = translation of the Phecode (for ICD-9 code, or protein name: if Phecode is the uniprot name, then Description can be the gene name, for example) Groupnum = the group number for plotting. If the code or analyte is part of a larger group (i.g. all lipid metabolites), then group these together and label by number here Group = identifier for the group (i.g. if Lipid metabolites are all group 1, then "Group" is the label, "Lipid Metabolite.") Color = R-palette color for plotting by group.


SNP_sex_interactions.r
Assess for Gene x sex interactions in Arivale data The PheWAS package returns a main effect beta coefficient only. For interaction analyses, I want the ability to look at ALL the coefficients. And specifically, I want to run FDR adjustment on the interaction term coefficients only. This code will run regression with SNP-sex interaction terms for every analyte. Chemistries need to be run separately from proteins and metabolites due to vendor effects. Results will be joined before FDR adjustment.


Boxplots_Arivale.r
Create boxplots of arivale analytes by genotype and decade This code uses ggplot to create a box plot of analyte values by genotype, across decades. SNP genotypes are color-coded, y-axis is log-transformed analyte value.
Also: creates box plots, faceted by sex.
Upload data file containing ALL the log transformed analyses, plus covariates & PCs info that was created in the SNP_sex_interactions_r.txt code

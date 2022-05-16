library(SKAT)
library(openxlsx)

# set working directory for dev and debug
bdir = "/Users/williamlee/SharpLab/Data/DVA_SKAT-O"
setwd(bdir)

# generate EUR SSD files
Generate_SSD_SetID('PLINK_EHv5_input/DVA_EUR_EHv5_SKATO.bed', 'PLINK_EHv5_input/DVA_EUR_EHv5_SKATO.bim', 'PLINK_EHv5_input/DVA_EUR_EHv5_SKATO.fam', 'SetID_files/DVA_EUR_EHv5.SetID', 'SSD_files/DVA_EUR_EHv5.SSD', 'SSD_files/DVA_EUR_EHv5.INFO')
# generate AFR SSD files
Generate_SSD_SetID('PLINK_EHv5_input/DVA_AFR_EHv5_SKATO.bed', 'PLINK_EHv5_input/DVA_AFR_EHv5_SKATO.bim', 'PLINK_EHv5_input/DVA_AFR_EHv5_SKATO.fam', 'SetID_files/DVA_AFR_EHv5.SetID', 'SSD_files/DVA_AFR_EHv5.SSD', 'SSD_files/DVA_AFR_EHv5.INFO')
# generate AMR SSD files
Generate_SSD_SetID('PLINK_EHv5_input/DVA_AMR_EHv5_SKATO.bed', 'PLINK_EHv5_input/DVA_AMR_EHv5_SKATO.bim', 'PLINK_EHv5_input/DVA_AMR_EHv5_SKATO.fam', 'SetID_files/DVA_AMR_EHv5.SetID', 'SSD_files/DVA_AMR_EHv5.SSD', 'SSD_files/DVA_AMR_EHv5.INFO')

pheno_df <- read.table(header=T, quote = "", sep="\t", file="REGENIE_pheno_Medium.tsv", stringsAsFactors=FALSE)[,2:3]

eur_fam <- Read_Plink_FAM_Cov('PLINK_EHv5_input/DVA_EUR_EHv5_SKATO.fam', 'REGENIE_covar.tsv', Is.binary=TRUE)
eur_fam$Phenotype <- NULL
eur_fam <- merge(x=eur_fam, y=pheno_df, by="IID", all.x=TRUE)

afr_fam <- Read_Plink_FAM_Cov('PLINK_EHv5_input/DVA_AFR_EHv5_SKATO.fam', 'REGENIE_covar.tsv', Is.binary=TRUE)
afr_fam$Phenotype <- NULL
afr_fam <- merge(x=afr_fam, y=pheno_df, by="IID", all.x=TRUE)

amr_fam <- Read_Plink_FAM_Cov('PLINK_EHv5_input/DVA_AMR_EHv5_SKATO.fam', 'REGENIE_covar.tsv', Is.binary=TRUE)
amr_fam$Phenotype <- NULL
amr_fam <- merge(x=amr_fam, y=pheno_df, by="IID", all.x=TRUE)

eur_obj <- SKAT_Null_Model(pheno ~ gender + PC1 + PC2 + PC3 + PC4 + PC5 + InsertSize, out_type="C", data=eur_fam)
afr_obj <- SKAT_Null_Model(pheno ~ gender + PC1 + PC2 + PC3 + PC4 + PC5 + InsertSize, out_type="C", data=afr_fam)
amr_obj <- SKAT_Null_Model(pheno ~ gender + PC1 + PC2 + PC3 + PC4 + PC5 + InsertSize, out_type="C", data=amr_fam)

SSD.INFO <- Open_SSD('SSD_files/DVA_EUR_EHv5.SSD', 'SSD_files/DVA_EUR_EHv5.INFO')
skato <- SKAT.SSD.All(SSD.INFO, eur_obj, method="SKATO")
Close_SSD()
eur_skato_df = skato[['results']]
eur_skato_df = eur_skato_df[order(eur_skato_df['P.value']),]
write.xlsx(eur_skato_df, file='out/EUR_EHv5_SKATO_out.xlsx')

SSD.INFO <- Open_SSD('SSD_files/DVA_AFR_EHv5.SSD', 'SSD_files/DVA_AFR_EHv5.INFO')
skato <- SKAT.SSD.All(SSD.INFO, afr_obj, method="SKATO")
Close_SSD()
afr_skato_df = skato[['results']]
afr_skato_df = afr_skato_df[order(afr_skato_df['P.value']),]
write.xlsx(afr_skato_df, file='out/AFR_EHv5_SKATO_out.xlsx')

SSD.INFO <- Open_SSD('SSD_files/DVA_AMR_EHv5.SSD', 'SSD_files/DVA_AMR_EHv5.INFO')
skato <- SKAT.SSD.All(SSD.INFO, amr_obj, method="SKATO")
Close_SSD()
amr_skato_df = skato[['results']]
amr_skato_df = amr_skato_df[order(amr_skato_df['P.value']),]
write.xlsx(amr_skato_df, file='out/AMR_EHv5_SKATO_out.xlsx')
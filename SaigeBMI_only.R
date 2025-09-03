### SaigeBMI_only as per new Regenie

library(dplyr)
head(OA_Cov_Regenie6Bmod)
dim(OA_Cov_Regenie6Bmod)  ### 487409     16

SaigeBMI_only<-OA_Cov_Regenie6Bmod%>%
                rename(Age = AgeatchronicPainQ,Sex = Gender)

dim(SaigeBMI_only)

write_csv(SaigeBMI_only, "SaigeBMI_only.csv")


#### Now, for SAIGE, we need to intersect Pain Phenotype together in the same dataset 
## with the covariates ### so we need to go back and retrieve the pheno file from Regenie ##

library(readr)
library(data.table)


samp = fread('template.txt', header = F, col.names = c('FID','IID','missing','sex')) %>%
  select(FID,IID)


dim(samp)  ### 487409      2

head(samp)  ### FID     IID

class(samp$FID)
class(samp$IID)


out6B = samp %>%
  left_join(OA_Hip_and_Knee_CaseControl3, by = 'FID')

dim(out6B)  ### 487409      3



out6B<-out6B%>%
  drop_na(FID,IID)

head(out6B)
dim(out6B)  ### 487409    3

class(out6B$FID)
class(out6B$IID)

out6B<-unique(out6B, by ="FID")

dim(out6B)  ### 487409      3

duplicated_entries<-out6B[out6B$FID %in% out6B$FID [duplicated(out6B$FID)] , ]

print(duplicated_entries)




dim(duplicated_entries)  ### 0 3

out6B<-out6B %>%
  filter(FID !=1794663)

dim(out6B)  ### 487408      3

write_csv(out6B, 'pain.phen6B.csv')

dim(pain_phen6B)  ### 487408      3
pain_phen6B

library(dplyr)


PhenoSaigeBMI_only<-left_join(pain_phen6B, SaigeBMI_only,by = "FID")

dim(PhenoSaigeBMI_only)  ### 487408     18

head(PhenoSaigeBMI_only)

PhenoSaigeBMI_only<-PhenoSaigeBMI_only[ , -c(4)]

dim(PhenoSaigeBMI_only)

PhenoSaigeBMI_only<-PhenoSaigeBMI_only%>%
                   rename(IID =IID.x)

write_csv(PhenoSaigeBMI_only,"PhenoSaigeBMI_only.csv")

write.table(PhenoSaigeBMI_only,"/rds/projects/s/sharmaoa-oapain/Saige_Stage1/SaigeBMI_only/PhenoSaigeBMI_only.tsv", sep = "\t", row.names = FALSE, quote = FALSE)










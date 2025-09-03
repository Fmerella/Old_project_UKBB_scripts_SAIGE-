ls()
load(Saige1_OAPain.rda)
setwd("/rds/projects/s/sharmaoa-oapain/Saige_Stage1")
load("/rds/projects/s/sharmaoa-oapain/Saige_Stage1/Saige1_OAPain.rda")
ls()
str(modglmm)
list_name[["modglmm"]]
modglmm[["coefficients"]]
modglmm[["linear.predictors"]]
         
library(readr)
library(data.table)

head(dragen_pvcf_coordinates.csv)
getwd()

vcf_coordinates<-read.csv("dragen_pvcf_coordinates.csv", header =TRUE)

library(dplyr)

vcf_files_b5103<-vcf_coordinates %>%
                            filter(filename =="ukb24310_c7_b5103_v1.vcf.gz")
                                                   
####, "ukb24310_c7_b5104_v1.vcf.gz", 
####               "ukb24310_c7_b5105_v1.vcf.gz", "ukb24310_c7_b5106_v1.vcf.gz", 
####                  "ukb24310_c7_b5107_v1.vcf.gz", "ukb24310_c7_b5108_v1.vcf.gz", 
###               "ukb24310_c7_b5109_v1.vcf.gz","ukb24310_c7_b511_v1.vcf.gz", 
###                "ukb24310_c7_b5110_v1.vcf.gz"))

write_csv(vcf_files_b5103,"vcf_files_b5103.csv")



vcf_files_b5104<-vcf_coordinates %>%
  filter(filename =="ukb24310_c7_b5104_v1.vcf.gz")

write_csv(vcf_files_b5104, "vcf_files_b5104.csv")


### Tomorrow we are going to focus on these codes based on the file obtained above which has been saved
## on Nexius #### and we are going to generate a csi index for SAIGE2 


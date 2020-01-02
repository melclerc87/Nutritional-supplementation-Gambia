rm(list = ls())

## LRTI Gambia supplementation analysis ##
## Nutrition 

# Load in data
library(data.table); library(vegan); library(fpc); library(gridExtra); library(ggplot2); library(Hmisc); library(ggdendro); library(tidyr); library(plyr); library(dplyr); library(grid); library(phyloseq); library(lme4);library(lmerTest); library(metagenomeSeq); library(EnhancedVolcano); library(VSURF); library(RVAideMemoire); library(multcompView); library(kableExtra); library(cowplot); library(pander); library(ggfortify); library(ResourceSelection); library(questionr); library(multcomp); library(lsmeans); library(afex); library(car); library(effects); library(nnet); library(MASS); library(xlsx)

load("U:/Datastore/CMVM/scs/groups/CIR-Bogaert_group/Melanie/09_project_Supplementation_LRTI_Gambia/analysis/ps_decontam.RData")

otu_table_gamb<-as(otu_table(ps_decontam$RA_decontam$sup_lrti_samples_final), "matrix")
otu_table_gamb_raw<-as(otu_table(ps_decontam$raw_decontam$sup_lrti_samples_final), "matrix")

meta_gamb<-data.frame(sample_data(ps_decontam$RA_decontam$sup_lrti_samples_final), check.names = F)%>%
  mutate_at(c("trt_grp", "cgroup", "viruspos", "viruscoinf", "rsv",  "stunting"), as.factor)%>%
  mutate(age_fac=cut(age, c(0,1,2,3,4,5), labels = c("1","2","3","4","5")), # each factor level contains the ages of that factor level and below
         shannon=as.numeric(vegan::diversity(t(otu_table_gamb), "shannon")),
         chao=estimateR(t(otu_table(ps_decontam$raw_decontam$sup_lrti_samples_final)))[2,],
         observed=colSums(ifelse(otu_table_gamb>0,1,0))) 

meta_gamb_d5<-meta_gamb%>%filter(timepoint=="D5")
meta_gamb_w6<-meta_gamb%>%filter(timepoint=="W6")
meta_gamb_d1<-meta_gamb%>%filter(timepoint=="D1")

meta_gamb_ext<-read.xlsx("U:/Datastore/CMVM/scs/groups/CIR-Bogaert_group/Melanie/09_project_Supplementation_LRTI_Gambia/metadata/Limited MMCT_LRTI data for Microbiome analysis with week 6 anthropometry.xlsx", sheetIndex = 4) # re-measurement values


# Effect of supplementation on stunting
table(meta_gamb_d5$trt_grp, meta_gamb_d5$stunting)

table(meta_gamb_ext$stunting, meta_gamb_ext$stunting3) # first argument = rows, second argument = columns
# over the 6-week period, two children recovered from stunting, one child became stunted

# Effect of supplementation on wasting
table(meta_gamb_ext$wasting, meta_gamb_ext$wasting3)
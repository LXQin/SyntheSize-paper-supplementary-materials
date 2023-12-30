library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggh4x)
library(RColorBrewer)
library(ggsci)
library(beanplot)
library(ggforce)
library(cowplot)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##################################
###### RESULTS ORGANIZATION ######
##################################

#### SKCM ####################### 
df_skcm <- skcm[skcm$pilot_size %in% c(20, 40, 60, 80, 100),]
df_skcm_AEhead = skcm_AEhead[skcm_AEhead$pilot_size %in% c(20, 60, 100) & skcm_AEhead$modelname != "MAF",]
df_skcm_epoch = skcm_epoch[skcm_epoch$epoch != "Epoch=300",]
df_skcm_norm = skcm_norm
df_skcm_batch = skcm_batch
rm(skcm)
rm(skcm_AEhead)
rm(skcm_epoch)
rm(skcm_norm)
rm(skcm_batch)

## SKCM check
skcm <- df_skcm
table(skcm$pilot_size, skcm$model)
save(skcm, file = "../FinalResAug2023/SKCM.RData")

## AEHead check
skcm_AEhead$modelname <- factor(skcm_AEhead$modelname, levels = c("glow","maf","realnvp"), labels =c("GLOW","MAF","RealNVP") )
skcm_AEhead <- rbind(skcm_AEhead, df_skcm_AEhead)
table(skcm_AEhead$pilot_size, skcm_AEhead$modelname, skcm_AEhead$offline)
save(skcm_AEhead, file = "../FinalResAug2023/SKCM_AEhead.RData")


## Batch check
skcm_batch$modelname = "MAF"
skcm_batch <- rbind(skcm_batch, df_skcm_batch)
table(skcm_batch$modelname, skcm_batch$pilot_size, skcm_batch$batch_frac)
save(skcm_batch, file = "../FinalResAug2023/SKCM_Batch.RData")

## Epoch check
skcm_epoch$modelname = "MAF"
skcm_epoch$epoch <- ifelse(skcm_epoch$epoch=="Early Stopping","EarlyStop","Fixed")
df_skcm_epoch$epoch <- ifelse(df_skcm_epoch$epoch=="EarlyStop","EarlyStop","Fixed")

skcm_epoch <- rbind(skcm_epoch, df_skcm_epoch)
table(skcm_epoch$modelname, skcm_epoch$pilot_size, skcm_epoch$epoch)
save(skcm_epoch, file = "../FinalResAug2023/SKCM_Epoch.RData")

## Normalization check
skcm_norm$modelname = "MAF"
skcm_norm <- rbind(skcm_norm, df_skcm_norm)
table(skcm_norm$modelname, skcm_norm$pilot_size, skcm_norm$normalization)
save(skcm_norm, file = "../FinalResAug2023/SKCM_Normalization.RData")


#### SKCM Positive #########
df_skcm_pos <- skcm_pos[skcm_pos$pilot_size %in% c(20, 40, 60, 80, 100),]
df_skcm_AEhead_pos = skcm_pos_AEhead[skcm_pos_AEhead$pilot_size %in% c(20, 60, 100) & skcm_pos_AEhead$modelname != "MAF",]
df_skcm_epoch_pos = skcm_pos_epoch[skcm_pos_epoch$epoch != "Epoch=300",]
df_skcm_norm_pos = skcm_pos_norm
df_skcm_batch_pos = skcm_pos_batch
rm(skcm_pos)
rm(skcm_pos_AEhead)
rm(skcm_pos_epoch)
rm(skcm_pos_norm)
rm(skcm_pos_batch)

## SKCMPos check
skcm_pos <- df_skcm_pos
table(skcm_pos$pilot_size, skcm_pos$model)
save(skcm_pos, file = "../FinalResAug2023/SKCMPositive_4.RData")
## AEHead check
skcm_pos_AEhead$modelname <- factor(skcm_pos_AEhead$modelname, levels = c("glow","maf","realnvp"), labels =c("GLOW","MAF","RealNVP") )
skcm_pos_AEhead <- rbind(skcm_pos_AEhead, df_skcm_AEhead_pos)
table(skcm_pos_AEhead$pilot_size, skcm_pos_AEhead$modelname, skcm_pos_AEhead$offline)
save(skcm_pos_AEhead, file = "../FinalResAug2023/SKCMPositive_4_AEhead.RData")


## Batch check
skcm_pos_batch$modelname = "MAF"
skcm_pos_batch <- rbind(skcm_pos_batch, df_skcm_batch_pos)
table(skcm_pos_batch$modelname, skcm_pos_batch$pilot_size, skcm_pos_batch$batch_frac)
save(skcm_pos_batch, file = "../FinalResAug2023/SKCMPositive_4_Batch.RData")

## Epoch check
skcm_pos_epoch$modelname = "MAF"
skcm_pos_epoch$epoch <- ifelse(skcm_pos_epoch$epoch=="Early Stopping","EarlyStop","Fixed")
df_skcm_epoch_pos$epoch <- ifelse(df_skcm_epoch_pos$epoch=="EarlyStop","EarlyStop","Fixed")

skcm_pos_epoch <- rbind(skcm_pos_epoch, df_skcm_epoch_pos)
table(skcm_pos_epoch$modelname, skcm_pos_epoch$pilot_size, skcm_pos_epoch$epoch)
save(skcm_pos_epoch, file = "../FinalResAug2023/SKCMPositive_4_Epoch.RData")

## Normalization check
skcm_pos_norm$modelname = "MAF"
skcm_pos_norm <- rbind(skcm_pos_norm, df_skcm_norm_pos)
table(skcm_pos_norm$modelname, skcm_pos_norm$pilot_size, skcm_pos_norm$normalization)
save(skcm_pos_norm, file = "../FinalResAug2023/SKCMPositive_4_Normalization.RData")


#### BRCA ###############
df_brca <- brca
df_brca_AEhead = brca_AEhead
df_brca_epoch = brca_epoch[brca_epoch$epoch != "Epoch=300",]
df_brca_norm = brca_norm
df_brca_batch = brca_batch
rm(brca)
rm(brca_AEhead)
rm(brca_epoch)
rm(brca_norm)
rm(brca_batch)

## brca check
brca$model <- factor(brca$model , levels = c("glow","maf","realnvp"), labels =c("GLOW","MAF","RealNVP") )
brca <- rbind(brca,df_brca)
table(brca$pilot_size, brca$model)
save(brca, file = "../FinalResAug2023/BRCA.RData")
## AEHead check
brca_AEhead$modelname <- factor(brca_AEhead$modelname, levels = c("glow","maf","realnvp"), labels =c("GLOW","MAF","RealNVP") )
brca_AEhead <- rbind(brca_AEhead, df_brca_AEhead)
table(brca_AEhead$pilot_size, brca_AEhead$modelname, brca_AEhead$offline)
save(brca_AEhead, file = "../FinalResAug2023/BRCA_AEhead.RData")


## Batch check
brca_batch$modelname = "MAF"
brca_batch <- rbind(brca_batch, df_brca_batch)
table(brca_batch$modelname, brca_batch$pilot_size, brca_batch$batch_frac)
save(brca_batch, file = "../FinalResAug2023/BRCA_Batch.RData")

## Epoch check
brca_epoch$modelname = "MAF"
brca_epoch$epoch <- ifelse(brca_epoch$epoch=="Early Stopping","EarlyStop","Fixed")
df_brca_epoch$epoch <- ifelse(df_brca_epoch$epoch=="EarlyStop","EarlyStop","Fixed")

brca_epoch <- rbind(brca_epoch, df_brca_epoch)
table(brca_epoch$modelname, brca_epoch$pilot_size, brca_epoch$epoch)
save(brca_epoch, file = "../FinalResAug2023/BRCA_Epoch.RData")

## Normalization check
brca_norm$modelname = "MAF"
brca_norm <- rbind(brca_norm, df_brca_norm)
table(brca_norm$modelname, brca_norm$pilot_size, brca_norm$normalization)
save(brca_norm, file = "../FinalResAug2023/BRCA_Normalization.RData")

#### BRCA Positive #########
df_brca_pos <- brca_pos
df_brca_AEhead_pos = brca_pos_AEhead
df_brca_epoch_pos = brca_pos_epoch[brca_pos_epoch$epoch != "Epoch=300",]
df_brca_norm_pos = brca_pos_norm
df_brca_batch_pos = brca_pos_batch
rm(brca_pos)
rm(brca_pos_AEhead)
rm(brca_pos_epoch)
rm(brca_pos_norm)
rm(brca_pos_batch)

## brca check
brca_pos$model <-  factor(brca_pos$model , levels = c("glow","maf","realnvp"), labels =c("GLOW","MAF","RealNVP") )
brca_pos <- rbind(brca_pos,df_brca_pos)
table(brca_pos$pilot_size, brca_pos$model)
save(brca_pos, file = "../FinalResAug2023/BRCAPositive_3.RData")

## AEHead check
brca_pos_AEhead$modelname <- factor(brca_pos_AEhead$modelname, levels = c("glow","maf","realnvp"), labels =c("GLOW","MAF","RealNVP") )
brca_pos_AEhead <- rbind(brca_pos_AEhead, df_brca_AEhead_pos)
table(brca_pos_AEhead$pilot_size, brca_pos_AEhead$modelname, brca_pos_AEhead$offline)
save(brca_pos_AEhead, file = "../FinalResAug2023/BRCAPositive_3_AEhead.RData")


## Batch check
brca_pos_batch$modelname = "MAF"
brca_pos_batch <- rbind(brca_pos_batch, df_brca_batch_pos)
table(brca_pos_batch$modelname, brca_pos_batch$pilot_size, brca_pos_batch$batch_frac)
save(brca_pos_batch, file = "../FinalResAug2023/BRCAPositive_3_Batch.RData")

## Epoch check
brca_pos_epoch$modelname = "MAF"
brca_pos_epoch$epoch <- ifelse(brca_pos_epoch$epoch=="Early Stopping","EarlyStop","Fixed")
df_brca_epoch_pos$epoch <- ifelse(df_brca_epoch_pos$epoch=="EarlyStop","EarlyStop","Fixed")

brca_pos_epoch <- rbind(brca_pos_epoch, df_brca_epoch_pos)
table(brca_pos_epoch$modelname, brca_pos_epoch$pilot_size, brca_pos_epoch$epoch)
save(brca_pos_epoch, file = "../FinalResAug2023/BRCAPositive_3_Epoch.RData")

## Normalization check
brca_pos_norm$modelname = "MAF"
brca_pos_norm <- rbind(brca_pos_norm, df_brca_norm_pos)
table(brca_pos_norm$modelname, brca_pos_norm$pilot_size, brca_pos_norm$normalization)
save(brca_pos_norm, file = "../FinalResAug2023/BRCAPositive_3_Normalization.RData")




#### SKCMLAML ###########
df_skcmlaml <- skcmlaml
df_skcmlaml_AEhead = skcmlaml_AEhead
df_skcmlaml_epoch = skcmlaml_epoch[skcmlaml_epoch$epoch != "Epoch=300",]
df_skcmlaml_norm = skcmlaml_norm
df_skcmlaml_batch = skcmlaml_batch
rm(skcmlaml)
rm(skcmlaml_AEhead)
rm(skcmlaml_epoch)
rm(skcmlaml_norm)
rm(skcmlaml_batch)

## skcmlaml check
skcmlaml$model <- "MAF"
skcmlaml <- rbind(skcmlaml,df_skcmlaml)
table(skcmlaml$pilot_size, skcmlaml$model)
save(skcmlaml, file = "../FinalResAug2023/SKCMLAML.RData")
## AEHead check
skcmlaml_AEhead$modelname <- "MAF"
skcmlaml_AEhead <- rbind(skcmlaml_AEhead, df_skcmlaml_AEhead)
table(skcmlaml_AEhead$pilot_size, skcmlaml_AEhead$modelname, skcmlaml_AEhead$offline)
save(skcmlaml_AEhead, file = "../FinalResAug2023/SKCMLAML_AEhead.RData")


## Batch check
skcmlaml_batch$modelname = "MAF"
skcmlaml_batch <- rbind(skcmlaml_batch, df_skcmlaml_batch)
table(skcmlaml_batch$modelname, skcmlaml_batch$pilot_size, skcmlaml_batch$batch_frac)
save(skcmlaml_batch, file = "../FinalResAug2023/SKCMLAML_Batch.RData")

## Epoch check
skcmlaml_epoch$modelname = "MAF"
skcmlaml_epoch$epoch <- ifelse(skcmlaml_epoch$epoch=="Early Stopping","EarlyStop","Fixed")
df_skcmlaml_epoch$epoch <- ifelse(df_skcmlaml_epoch$epoch=="EarlyStop","EarlyStop","Fixed")

skcmlaml_epoch <- rbind(skcmlaml_epoch, df_skcmlaml_epoch)
table(skcmlaml_epoch$modelname, skcmlaml_epoch$pilot_size, skcmlaml_epoch$epoch)
save(skcmlaml_epoch, file = "../FinalResAug2023/SKCMLAML_Epoch.RData")

## Normalization WATINTING
skcmlaml_norm$modelname = "MAF"
skcmlaml_norm <- rbind(skcmlaml_norm, df_skcmlaml_norm)
table(skcmlaml_norm$modelname, skcmlaml_norm$pilot_size, skcmlaml_norm$normalization)
save(skcmlaml_norm, file = "../FinalResAug2023/SKCMLAML_Normalization.RData")

#### SKCMLAML Positive #########
df_skcmlaml_pos <- skcmlaml_pos
df_skcmlaml_AEhead_pos = skcmlaml_pos_AEhead
df_skcmlaml_epoch_pos = skcmlaml_pos_epoch[skcmlaml_pos_epoch$epoch != "Epoch=300",]
df_skcmlaml_norm_pos = skcmlaml_pos_norm
df_skcmlaml_batch_pos = skcmlaml_pos_batch
rm(skcmlaml_pos)
rm(skcmlaml_pos_AEhead)
rm(skcmlaml_pos_epoch)
rm(skcmlaml_pos_norm)
rm(skcmlaml_pos_batch)

## skcmlaml check
skcmlaml_pos$model <- "MAF"
skcmlaml_pos <- rbind(df_skcmlaml_pos,skcmlaml_pos)
table(skcmlaml_pos$pilot_size, skcmlaml_pos$model)
save(skcmlaml_pos, file = "../FinalResAug2023/SKCMLAMLPositive_3.RData")
## AEHead check
skcmlaml_pos_AEhead$modelname <- "MAF"
skcmlaml_pos_AEhead <- rbind(skcmlaml_pos_AEhead, df_skcmlaml_AEhead_pos)
table(skcmlaml_pos_AEhead$pilot_size, skcmlaml_pos_AEhead$modelname, skcmlaml_pos_AEhead$offline)
save(skcmlaml_pos_AEhead, file = "../FinalResAug2023/SKCMLAMLPositive_3_AEhead.RData")


## Batch check
skcmlaml_pos_batch$modelname = "MAF"
skcmlaml_pos_batch <- rbind(skcmlaml_pos_batch, df_skcmlaml_batch_pos)
table(skcmlaml_pos_batch$modelname, skcmlaml_pos_batch$pilot_size, skcmlaml_pos_batch$batch_frac)
save(skcmlaml_pos_batch, file = "../FinalResAug2023/SKCMLAMLlPositive_3_Batch.RData")

## Epoch check
skcmlaml_pos_epoch$modelname = "MAF"
skcmlaml_pos_epoch$epoch <- ifelse(skcmlaml_pos_epoch$epoch=="Early Stopping","EarlyStop","Fixed")
df_skcmlaml_epoch_pos$epoch <- ifelse(df_skcmlaml_epoch_pos$epoch=="EarlyStop","EarlyStop","Fixed")

skcmlaml_pos_epoch <- rbind(skcmlaml_pos_epoch, df_skcmlaml_epoch_pos)
table(skcmlaml_pos_epoch$modelname, skcmlaml_pos_epoch$pilot_size, skcmlaml_pos_epoch$epoch)
save(skcmlaml_pos_epoch, file = "../FinalResAug2023/SKCMLAMLPositive_3_Epoch.RData")

## Normalization check
skcmlaml_pos_norm$modelname = "MAF"
skcmlaml_pos_norm <- rbind(skcmlaml_pos_norm, df_skcmlaml_norm_pos)
table(skcmlaml_pos_norm$modelname, skcmlaml_pos_norm$pilot_size, skcmlaml_pos_norm$normalization)
save(skcmlaml_pos_norm, file = "../FinalResAug2023/SKCMLAMLPositive_3_Normalization.RData")




#### LAML Positive ####
laml_pos_AEhead$epoch <- "EarlyStop"
laml_pos_opt <- rbind(laml_pos_opt, laml_pos_AEhead)
table(laml_pos_opt$pilot_size, laml_pos_opt$modelname, laml_pos_opt$epoch)
save(laml_pos_opt, file = "../FinalResAug2023/LAMLPositive_2_Optimal.RData")

#### PRAD Positive ####
prad_pos_AEhead$epoch <- "EarlyStop"
prad_pos_opt <- rbind(prad_pos_opt, prad_pos_AEhead)
table(prad_pos_opt$pilot_size, prad_pos_opt$modelname, prad_pos_opt$epoch)
save(prad_pos_opt, file = "../FinalResAug2023/PRADPositive_2_Optimal.RData")

#### BRCAPRAD Positive ####
brcaprad_pos_AEhead$normalization <- "TC"
brcaprad_pos_opt <- rbind(brcaprad_pos_opt, brcaprad_pos_AEhead)
table(brcaprad_pos_opt$pilot_size, brcaprad_pos_opt$modelname, brcaprad_pos_opt$normalization)
save(brcaprad_pos_opt, file = "../FinalResAug2023/BRCAPRADPositive_3_Optimal.RData")


#### RNABRCAPos ####
# RNABRCA check
rnabrca_pos <- rbind(rnabrca_pos, RNABRCA_Pos)
table(rnabrca_pos$pilot_size, rnabrca_pos$model)
rnabrca_pos$modelname <- rnabrca_pos$model
rnabrca_pos$modelname[rnabrca_pos$model == "maf"] <- "MAF"
rnabrca_pos$modelname[rnabrca_pos$model == "epochES_WGANGP"] <- "WGANGP"
save(rnabrca_pos, file = "../FinalResAug2023/RNABRCAPositive_5-2.RData")
rm(rnabrca_pos)
rm(RNABRCA_Pos)
ggplot(rnabrca_pos)+geom_boxplot(aes(x = pilot_size, y=ks_zero, color = modelname))

# Gaussianhead check
int_func <- function(x){
  if(grepl( "Gaussianhead", strsplit(x, "_")[[1]][1],fixed = TRUE)){return("Gaussianhead")}
  else{return("Original")}}
RNABRCA_Pos_GH$offline <- sapply(RNABRCA_Pos_GH$model, int_func)
RNABRCA_Pos_GH$modelname <- "MAF"
rnabrca_pos_Gaussianhead <- rbind(rnabrca_pos_Gaussianhead, RNABRCA_Pos_GH)
table(rnabrca_pos_Gaussianhead$pilot_size, rnabrca_pos_Gaussianhead$modelname, rnabrca_pos_Gaussianhead$offline)
ggplot(rnabrca_pos_Gaussianhead)+geom_boxplot(aes(x = pilot_size, y=ks_zero, color = offline))+
  facet_wrap(vars(modelname))
save(rnabrca_pos_Gaussianhead, file = "../FinalResAug2023/RNABRCAPositive_5-2_Gaussianhead.RData")
rm(rnabrca_pos_Gaussianhead)
rm(RNABRCA_Pos_GH)

# Normalization check
table(RNABRCAPos_norm$pilot_size, RNABRCAPos_norm$normalization, RNABRCAPos_norm$modelname)
RNABRCAPos_norm$modelname <- "MAF"
rnabrca_pos_norm <- RNABRCAPos_norm
save(rnabrca_pos_norm, file = "../FinalResAug2023/RNABRCAPositive_5-2_Normalization.RData")
ggplot(rnabrca_pos_norm)+geom_boxplot(aes(x = pilot_size, y=ARI, color = normalization))
rm(rnabrca_pos_norm)
rm(RNABRCAPos_norm)

#### RNAPRADPos ####

# RNAPRAD check
rnaprad_pos <- rbind(rnaprad_pos, RNAPRAD_Pos)
table(rnaprad_pos$pilot_size, rnaprad_pos$model)
rnaprad_pos$modelname <- rnaprad_pos$model
rnaprad_pos$modelname[rnaprad_pos$model == "maf"] <- "MAF"
rnaprad_pos$modelname[rnaprad_pos$model == "epochES_WGANGP"] <- "WGANGP"
ggplot(rnaprad_pos)+geom_boxplot(aes(x = pilot_size, y=ks_sd, color = modelname))
save(rnaprad_pos, file = "../FinalResAug2023/RNAPRADPositive_5-2.RData")
rm(rnaprad_pos)
rm(RNAPRAD_Pos)

# Gaussianhead check
int_func <- function(x){
  if(grepl( "Gaussianhead", strsplit(x, "_")[[1]][1],fixed = TRUE)){return("Gaussianhead")}
  else{return("Original")}}
RNAPRAD_Pos_GH$offline <- sapply(RNAPRAD_Pos_GH$model, int_func)
RNAPRAD_Pos_GH$modelname <- "MAF"
rnaprad_pos_Gaussianhead <- rbind(rnaprad_pos_Gaussianhead, RNAPRAD_Pos_GH)
table(rnaprad_pos_Gaussianhead$pilot_size, rnaprad_pos_Gaussianhead$modelname, rnaprad_pos_Gaussianhead$offline)
ggplot(rnaprad_pos_Gaussianhead)+geom_boxplot(aes(x = pilot_size, y=ks_zero, color = offline))+
  facet_wrap(vars(modelname))
save(rnaprad_pos_Gaussianhead, file = "../FinalResAug2023/RNAPRADAPositive_5-2_Gaussianhead.RData")
rm(rnaprad_pos_Gaussianhead)
rm(RNAPRAD_Pos_GH)

# Normalization check
table(RNAPRADPos_norm$pilot_size, RNAPRADPos_norm$normalization, RNAPRADPos_norm$modelname)
RNAPRADPos_norm$modelname <- "MAF"
rnaprad_pos_norm <- RNAPRADPos_norm
save(rnaprad_pos_norm, file = "../FinalResAug2023/RNAPRADPositive_5-2_Normalization.RData")
ggplot(rnaprad_pos_norm)+geom_boxplot(aes(x = pilot_size, y=ks_zero, color = normalization))
rm(rnaprad_pos_norm)
rm(RNAPRADPos_norm)



#### RNABRCAPRADPos ####

# RNABRCAPRAD check
rnabrcaprad_pos <- rbind(rnabrcaprad_pos, RNABRCAPRAD_pos)
table(rnabrcaprad_pos$pilot_size, rnabrcaprad_pos$model)
rnabrcaprad_pos$modelname <- rnabrcaprad_pos$model
rnabrcaprad_pos$modelname[rnabrcaprad_pos$model == "maf"] <- "MAF"
ggplot(rnabrcaprad_pos)+geom_boxplot(aes(x = pilot_size, y=ks_sd, color = modelname))
save(rnabrcaprad_pos, file = "../FinalResAug2023/RNABRCAPRADPositive_5-2.RData")
rm(rnabrcaprad_pos)
rm(RNABRCAPRAD_pos)

# Gaussianhead check
int_func <- function(x){
  if(grepl( "Gaussianhead", strsplit(x, "_")[[1]][1],fixed = TRUE)){return("Gaussianhead")}
  else{return("Original")}}
RNABRCAPRAD_GH$offline <- sapply(RNABRCAPRAD_GH$model, int_func)
RNABRCAPRAD_GH$modelname <- "MAF"
rnabrcaprad_pos_Gaussianhead$offline <- sapply(rnabrcaprad_pos_Gaussianhead$model, int_func)
rnabrcaprad_pos_Gaussianhead$modelname <- "CVAE1-100"
rnabrcaprad_pos_Gaussianhead <- rbind(rnabrcaprad_pos_Gaussianhead, RNABRCAPRAD_GH)
table(rnabrcaprad_pos_Gaussianhead$pilot_size, rnabrcaprad_pos_Gaussianhead$modelname, rnabrcaprad_pos_Gaussianhead$offline)
ggplot(rnabrcaprad_pos_Gaussianhead)+geom_boxplot(aes(x = pilot_size, y=ARI, color = offline))+
  facet_wrap(vars(modelname))
save(rnabrcaprad_pos_Gaussianhead, file = "../FinalResAug2023/RNABRCAPRADPositive_5-2_Gaussianhead.RData")
rm(rnabrcaprad_pos_Gaussianhead)
rm(RNABRCAPRAD_GH)

# Normalization check
table(RNABRCAPRAD_norm$pilot_size, RNABRCAPRAD_norm$normalization, RNABRCAPRAD_norm$modelname)
RNABRCAPRAD_norm$modelname <- "MAF"
rnabrcaprad_pos_norm <- RNABRCAPRAD_norm
save(rnabrcaprad_pos_norm, file = "../FinalResAug2023/RNABRCAPRADPositive_5-2_Normalization.RData")
ggplot(rnabrcaprad_pos_norm)+geom_boxplot(aes(x = pilot_size, y=ks_sd, color = normalization))
rm(rnabrcaprad_pos_norm)
rm(RNABRCAPRAD_norm)



#### BRCA Transfer ####
brca_transfer_my <- brca_transfer
brca_transfer$modelname <- "MAF"
brca_transfer <- rbind(brca_transfer, brca_transfer_my)
table(brca_transfer$modelname, brca_transfer$transfer, brca_transfer$pilot_size)
save(brca_transfer, file = "../FinalResAug2023/BRCA_Transfer.RData")
rm(brca_transfer, brca_transfer_my)
#### SKCM Transfer #### 
skcm_transfer_my <- skcm_transfer
skcm_transfer$modelname <- "MAF"
skcm_transfer <- rbind(skcm_transfer, skcm_transfer_my)
table(skcm_transfer$modelname, skcm_transfer$transfer, skcm_transfer$pilot_size)
save(skcm_transfer, file = "../FinalResAug2023/SKCM_Transfer.RData")
rm(skcm_transfer, skcm_transfer_my)
##################################
######## VISUALIZATION ###########
##################################

#### Full evaluation vis ####
vis_all <- function(dataname, nickname, cancer_type = c(1, 2)){
  # This function visualization the results of DGM evaluation
  # dataname = "SKCMPositive_4", nickname = "skcm_pos", cancer_type = 1
  
  
  # Step1: load data ####
  default <- "~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Interns/2021 Summer/Generative Models/DGMmicroRNA/FinalResAug2023/"
  load(paste(default, dataname, ".RData", sep = ""))
  load(paste(default, dataname, "_AEhead.RData", sep = ""))
  load(paste(default, dataname, "_Normalization.RData", sep = ""))
  load(paste(default, dataname, "_Batch.RData", sep = ""))
  load(paste(default, dataname, "_Epoch.RData", sep = ""))
  overall <- get(nickname)
  offline <- get(paste(nickname, "_AEhead", sep = ""))
  norm <- get(paste(nickname, "_norm", sep = ""))
  batch <-  get(paste(nickname, "_batch", sep = ""))
  epoch <-  get(paste(nickname, "_epoch", sep = ""))
  
  # Step2: generate dataframes for plotting  ####
  if(cancer_type == 1){
    # Evaluation criteria: 1-Prop.successful.genes, 2-CCC.pos.control, 3-cARI, 4-MAD
    mycolors <- c(brewer.pal(n = 8, name = "Blues")[c(2,4,6)], brewer.pal(n = 8, name = "Oranges")[c(2,4,6)], brewer.pal(n = 8, name = "Greens")[c(2,4,6)])
    plot_width <- c(10,10,16,7,7,7)
    plot_height <- c(6,6,6,6,6,6)
    df_proc <- function(res){
      res$succ_prop <- 1 - res$fail_prop
      res$pilot_size <- factor(res$pilot_size, levels = as.character(sort(as.numeric(as.character(unique(res$pilot_size))))))
      if("modelname" %in% colnames(res)){model <- res[,"modelname"]}
      else{model <- res[,"model"]}
      if(length(unique(model)) <= 3){lvs <- c("VAE1-10","MAF","WGANGP")}
      else{lvs <- c("VAE1-1", "VAE1-5", "VAE1-10","RealNVP","GLOW","MAF","GAN", "WGAN", "WGANGP")}
      res$model <- factor(model, levels = lvs)
      res$modeltype <- rep(NA, nrow(res))
      for (i in 1:nrow(res)) {
        if(substr(res$model[i], 1,3)=="VAE"){res$modeltype[i] <- "VAEs"}
        else if(substr(res$model[i], 1,3) %in% c("GAN", "WGA")){res$modeltype[i] <- "GANs"}
        else{res$modeltype[i] <- "FLOWs"}
      }
      res$modeltype <- factor(res$modeltype, levels = c("VAEs","FLOWs","GANs"))
      res$ccc_pos[res$ccc_pos<0 | is.na(res$ccc_pos)] <- 0
      res$ARI[res$ARI<0 | is.na(res$ARI)] <- 0
      res <- tidyr::pivot_longer(res, cols = c("ARI", "ccc_pos", "succ_prop"), names_to = "Metric", values_to = "Rate")
      res$Metric <- factor(res$Metric, levels = c("succ_prop","ARI","ccc_pos"),
                           labels = c("1 - Pct(0-genes)", "cARI", "CCC of PCC"))
      if("normalization" %in% colnames(res)){res$normalization <- factor(res$normalization, levels=c('None','TC','TMM','UQ'))}
      return(res)
    }
    heatmap_proc <- function(res){
      if("modelname" %in% colnames(res)){model <- res[,"modelname"]}
      else{model <- res[,"model"]}
      if(length(unique(model)) <= 4){lvs <- c("VAE1-10","MAF","WGANGP")}
      else{lvs <- c("VAE1-1", "VAE1-5", "VAE1-10","RealNVP","GLOW","MAF","GAN", "WGAN", "WGANGP")}
      res$model <- factor(model, levels = rev(lvs))
      res$pilot_size <- factor(res$pilot_size, levels = as.character(sort(as.numeric(as.character(unique(res$pilot_size))))))
      res$modeltype <- rep(NA, nrow(res))
      for (i in 1:nrow(res)) {
        if(substr(res$model[i], 1,3)=="VAE"){res$modeltype[i] <- "VAEs"}
        else if(substr(res$model[i], 1,3) %in% c("GAN", "WGA")){res$modeltype[i] <- "GANs"}
        else{res$modeltype[i] <- "FLOWs"}
      }
      res$modeltype <- factor(res$modeltype, levels = c("VAEs","FLOWs","GANs"))
      res <- res %>% group_by(pilot_size, model, modeltype) %>%
        summarise(
          Mean = mean(ks_mean, na.rm=TRUE),
          Sd = mean(ks_sd, na.rm=TRUE),
          Sparsity = mean(ks_zero, na.rm=TRUE)
        )
      
      res <- tidyr::pivot_longer(res, cols = c("Mean", "Sd", "Sparsity"), names_to = "Metric", values_to = "Rate")
      res$Metric <- factor(res$Metric, levels = (c("Mean", "Sd", "Sparsity")), labels = c("MAD(Mean)","MAD(SD)","MAD(Sparsity)"))
      res$text <- round(res$Rate, 3)
      res$text[res$Rate > 6 | is.na(res$Rate)] <- ">6"
      res$Rate[res$Rate > 6 | is.na(res$Rate)] <- 6
      res$modelorder <-  rep(NA, nrow(res))
      for (i in 1:nrow(res)) {
        if(res$model[i] %in% c("VAE1-1","GAN","RealNVP")){res$modelorder[i] <- "VAE1-1/RealNVP/GAN"}
        else if(res$model[i] %in% c("VAE1-5","WGAN","GLOW")){res$modelorder[i] <- "VAE1-5/GLOW/WGAN"}
        else{res$modelorder[i] <- "VAE1-10/MAF/WGANGP"}
      }
      res$modelorder <- factor(res$modelorder, levels = c("VAE1-10/MAF/WGANGP", "VAE1-5/GLOW/WGAN", "VAE1-1/RealNVP/GAN"))
      if("normalization" %in% colnames(res)){res$normalization <- factor(res$normalization, levels=c('None','TC','TMM','UQ'))}
      return(res)
    }
  }
  else{
    # Evaluation criteria: 1-Prop.successful.genes, 2-CCC.pos.control, 3-ARI, 4-MAD, 5-ccc.log10pvalue, 6-ccc.log2FC
    plot_width <- c(6,4,8,6,6,6)
    plot_height <- c(9.5, 9.5, 9.5, 9.5, 9.5, 9.5)
    mycolors <- c(brewer.pal(n = 8, name = "Blues")[c(2,4,6)],brewer.pal(n = 8, name = "Oranges")[6])
    df_proc <- function(res){
      res$succ_prop <- 1 - res$fail_prop
      res$pilot_size <- factor(res$pilot_size, levels = as.character(sort(as.numeric(as.character(unique(res$pilot_size))))))
      if("modelname" %in% colnames(res)){model <- res[,"modelname"]}
      else{model <- res[,"model"]}
      if(length(unique(model)) <= 2){lvs <- c("CVAE1-10","MAF")}
      else{lvs <- c("CVAE1-1", "CVAE1-5", "CVAE1-10","MAF")}
      res$model <- factor(model, levels = lvs)
      res$modeltype <- rep(NA, nrow(res))
      for (i in 1:nrow(res)) {
        if(substr(res$model[i], 1,3)=="CVA"){res$modeltype[i] <- "VAEs"}
        else{res$modeltype[i] <- "FLOWs"}
      }
      res$modeltype <- factor(res$modeltype, levels = c("VAEs","FLOWs"))
      res$ccc_pos[res$ccc_pos<0 | is.na(res$ccc_pos)] <- 0
      res$ARI[res$ARI<0 | is.na(res$ARI)] <- 0
      res$ccc_log10pvalue[res$ccc_log10pvalue < 0 | is.na(res$ccc_log10pvalue)] <- 0
      res$ccc_log2FC[res$ccc_log2FC < 0 | is.na(res$ccc_log2FC)] <- 0
      if("normalization" %in% colnames(res)){res$normalization <- factor(res$normalization, levels=c('None','TC','TMM','UQ','DESeq'))}
      res <- tidyr::pivot_longer(res, cols = c("ccc_log10pvalue", "ccc_log2FC","ARI", "succ_prop","ccc_pos"), names_to = "Metric", values_to = "Rate")
      res$Metric <- factor(res$Metric, levels = c("succ_prop","ARI","ccc_pos","ccc_log10pvalue", "ccc_log2FC"),
                           labels = c("1 - Pct(0-genes)", "ARI","CCC of PCC", "CCC of -log10(pvalue)", "CCC of log2FC"))
      return(res)
    }
    
    heatmap_proc <- function(res){
      if("modelname" %in% colnames(res)){model <- res[,"modelname"]}
      else{model <- res[,"model"]}
      if(length(unique(model)) <= 2){lvs <- c("CVAE1-10","MAF")}
      else{lvs <- c("CVAE1-1", "CVAE1-5", "CVAE1-10","MAF")}
      res$model <- factor(model, levels = lvs)
      res$pilot_size <- factor(res$pilot_size, levels = as.character(sort(as.numeric(as.character(unique(res$pilot_size))))))
      res$modeltype <- rep(NA, nrow(res))
      for (i in 1:nrow(res)) {
        if(substr(res$model[i], 1,3)=="CVA"){res$modeltype[i] <- "VAEs"}
        else{res$modeltype[i] <- "FLOWs"}
      }
      res$modeltype <- factor(res$modeltype, levels = c("VAEs","FLOWs"))
      res <- res %>% group_by(pilot_size, model, modeltype) %>%
        summarise(
          Mean = mean(ks_mean,na.rm=TRUE),
          Sd = mean(ks_sd,na.rm=TRUE),
          Sparsity = mean(ks_zero,na.rm=TRUE)
        )
      
      res <- tidyr::pivot_longer(res, cols = c("Mean", "Sd", "Sparsity"), names_to = "Metric", values_to = "Rate")
      res$Metric <- factor(res$Metric, levels = (c("Mean", "Sd", "Sparsity")), labels = c("MAD(Mean)","MAD(SD)","MAD(Sparsity)"))
      res$text <- round(res$Rate, 3)
      res$text[res$Rate > 6 | is.na(res$Rate)] <- ">6"
      res$Rate[res$Rate > 6 | is.na(res$Rate)] <- 6
      res$modelorder <- factor(res$model, levels=rev(levels(res$model)))
      res$modeltype <- "All models"
      if("normalization" %in% colnames(res)){res$normalization <- factor(res$normalization, levels=c('None','TC','TMM','UQ','DESeq'))}
      return(res)
    }
  }
  # Step 3: make graphs  ####
  gg_default <- function(p){
    pp <- p +
      theme(legend.position = "bottom")+
      # guides(fill=guide_legend(nrow = 1,byrow = TRUE))+
      theme(strip.text.x = element_text(face = "bold", size = 12),
            strip.text.y = element_text(face = "bold", size = 10),
            legend.title=element_text(size = 10), 
            legend.text=element_text(size = 9),
            axis.text = element_text(size = 10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank())
    return(pp)
  }
  #### p11 overall boxplot ####
  res <- df_proc(overall)
  p11 <- ggplot(data = res)+ 
    geom_violin(aes(x = pilot_size, y = Rate, fill = model), position = position_dodge(0.6), width=1)+
    # geom_boxplot(aes(x = pilot_size, y = Rate, fill = model), show.legend = T, width = 0.15,position = position_dodge(0.9))+
    facet_grid(row = vars(Metric), col = vars(modeltype))+
    theme_bw()+
    labs(x = "Pilot size", y = "", switch = "y", fill = "Models")+
    scale_fill_manual(values = mycolors)+
    theme(panel.border = element_rect(colour = "black", fill = NA))+
    ylim(0,1)+
    guides(fill=guide_legend(nrow = 1,byrow = TRUE))
  pdf(file = paste(default, dataname, ".pdf", sep = ""), width = plot_width[1], height = plot_height[1])
  print(gg_default(p11))
  dev.off()
  #### p12 overall heatmap ####
  res <- heatmap_proc(overall)
  p12 <- ggplot(res, aes(x = pilot_size , y = modelorder, fill = Rate)) + 
    geom_tile(show.legend = T)+
    labs(y = "",x = "Pilot size", fill = "MAD")+
    facet_nested(Metric~modeltype)+
    geom_text(aes(pilot_size, modelorder, label = text), color = "black", size = 4) +
    theme_minimal()+
    # scale_y_discrete(position = "left", labels=c( "Mean" = "MAD(Mean)", "Sd" = "MAD(SD)", "Sparsity" = "MAD(Sparsity)"))+
    scale_fill_gradientn(colours = terrain.colors(10))+
    theme(axis.text.y = element_text(angle = 0, hjust = 1, colour = "black", size = 8))+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank())
  pdf(file = paste(default, dataname, "_MAD.pdf", sep = ""), width = plot_width[2], height = plot_height[2])
  print(gg_default(p12))
  dev.off()
  
  #### p21 offline violine plot####
  res <- df_proc(offline)
  p21 <- ggplot(data = res, aes(x = pilot_size, y = Rate))+
    geom_violin(aes(color = offline), show.legend = T, position = position_dodge(0.9)) +
    geom_boxplot(aes(color = offline), width = 0.15, position = position_dodge(0.9)) +
    facet_grid(cols=vars(model), rows= vars(Metric), scale = "free_y")+
    theme_bw()+
    labs(x="Pilot size", y="", color = "Offline")+
    scale_y_continuous(position = "left")+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.line = element_blank())+
    scale_color_npg()+
    ylim(0,1)
  pdf(file = paste(default, dataname, "_Offline.pdf", sep = ""), width = plot_width[3], height = plot_height[3])
  print(gg_default(p21))
  dev.off()
  
  
  #### p3 batch sina plot ####
  res <- df_proc(batch)
  res$batch_size <- paste(res$batch_frac, "of pilot size")
  p3 <- ggplot(res[res$batch_frac!= "50%",])+ 
    # geom_sina(aes(x = pilot_size, y = Rate, color = batch_frac), show.legend = T)+
    geom_violin(aes(x = pilot_size, y = Rate, color = batch_size), show.legend = T, position = position_dodge(0.9)) +
    geom_boxplot(aes(x = pilot_size, y = Rate, color = batch_size), width = 0.15, position = position_dodge(0.9)) +
    facet_nested(Metric~model)+
    theme_bw()+
    labs(x="Pilot size", y="", color = "Batch size")+
    scale_y_continuous(position = "left")+
    theme(panel.border = element_rect(colour = "black", fill = NA),
          axis.line = element_blank())+
    scale_color_npg()+
    ylim(0,1)
  pdf(file = paste(default, dataname, "_Batch.pdf", sep = ""), width = plot_width[4], height = plot_height[4])
  print(gg_default(p3))
  dev.off()
  
  #### p4 epoch line plot ####
  # if(cancer_type == 1){lvs <- c("VAE1-10","MAF","WGANGP")}
  # else{lvs <- c("CVAE1-10","MAF")}
  # epoch$model <- factor(epoch$modelname, levels = lvs)
  # res <- epoch %>% group_by(model, pilot_size, epoch) %>%
  #   summarise(
  #     sd_ccc_pos = sd(ccc_pos,na.rm = T),
  #     mean_ccc_pos = mean(ccc_pos,na.rm = T),
  #     sd_ARI = sd(ARI,na.rm = T),
  #     mean_ARI = mean(ARI,na.rm = T),
  #   )
  # 
  # res_sd <- tidyr::pivot_longer(res, cols = c("sd_ccc_pos", "sd_ARI"), names_to = "Metric", values_to = "Sd")
  # res_sd <- data.frame(model = res_sd$model,pilot_size = res_sd$pilot_size, epoch = res_sd$epoch,
  #                      Metric = substr(res_sd$Metric, 4, 10), Sd = res_sd$Sd)
  # res_mean <- tidyr::pivot_longer(res, cols = c("mean_ccc_pos", "mean_ARI"), names_to = "Metric", values_to = "Mean")
  # res_mean <- data.frame(model = res_mean$model,pilot_size = res_mean$pilot_size, epoch = res_mean$epoch,
  #                        Metric = substr(res_mean$Metric, 6, 15), Mean = res_mean$Mean)
  # res <- inner_join(res_sd, res_mean, by = c("model","epoch","pilot_size","Metric"))
  # res$Metric <- factor(res$Metric, levels = c("ccc_pos","ARI"),
  #                      labels = c("CCC of PCC", "cARI"))
  # p4 <- ggplot(res, aes(x = pilot_size, y = Mean, group = epoch, color = epoch))+ 
  #   geom_errorbar(aes(ymin = Mean-Sd, ymax = Mean+Sd), width = 4) +
  #   geom_line() + 
  #   geom_point()+
  #   theme_bw()+
  #   facet_nested(Metric~model,scales = "free_y")+
  #   labs(x="Pilot size", y="", color = "Epoch")+
  #   scale_y_continuous(position = "left")+
  #   theme(panel.border = element_rect(colour = "black", fill = NA),
  #         axis.line = element_blank())+
  #   scale_color_npg()
  res <- df_proc(epoch)
  p4 <- ggplot(res)+ 
    geom_violin(aes(x = pilot_size, y = Rate, color = epoch), show.legend = T, position = position_dodge(0.9)) +
    geom_boxplot(aes(x = pilot_size, y = Rate, color = epoch), width = 0.15, position = position_dodge(0.9)) +
    facet_nested(Metric~model)+
    theme_bw()+
    labs(x="Pilot size", y="", color = "Epoch")+
    scale_y_continuous(position = "left")+
    theme(panel.border = element_rect(colour = "black", fill = NA),
          axis.line = element_blank())+
    scale_color_npg()+
    ylim(0,1)
  pdf(file = paste(default, dataname, "_Epoch.pdf", sep = ""), width = plot_width[5], height = plot_height[5])
  print(gg_default(p4))
  dev.off()
  
  #### p5 normalization violin plot ####
  res <- df_proc(norm)
  p5 <- ggplot(data = res[res$normalization!="DESeq",])+
    geom_violin(aes(x = pilot_size, y = Rate, fill = normalization ), position = position_dodge(0.9))+
    facet_grid(rows = vars(Metric), cols = vars(model),scales = "free_y")+
    theme_bw()+
    labs(x="Pilot size", y="", fill="Normalization")+
    scale_y_continuous(position = "left")+
    theme(panel.border = element_rect(colour = "black", fill = NA),
          axis.line = element_blank())+
    scale_fill_npg()+
    ylim(0,1)
  pdf(file = paste(default, dataname, "_Normalization.pdf", sep = ""), width = plot_width[6], height = plot_height[6])
  print(gg_default(p5))
  dev.off()
  
}

vis_all("SKCM", "skcm", 1)
vis_all("SKCMPositive_4", "skcm_pos", 1)
vis_all("BRCA", "brca", 1)
vis_all("BRCAPositive_3", "brca_pos", 1)
vis_all("SKCMLAML", "skcmlaml", 2)
vis_all("SKCMLAMLPositive_3", "skcmlaml_pos", 2)

#### Two-group UMAP visualization #### 

# run ggdefault in vis_all(), UMAP_eval in Evaluation.rmd first
default <- "~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Interns/2021 Summer/Generative Models/DGMmicroRNA/FinalResAug2023/"
selected_pilot <- c("20", "60", "100")
selected_model <- c("CVAE1-10")
g1 <- "SKCM"
g2 <- "LAML"
dataname <- "SKCMLAMLPositive_3"
real <- read.csv(paste("../RealData/", dataname, ".csv", sep = ""), header = T)
n1 <- sum(real$groups==g1)
n2 <- sum(real$groups==g2)
real <- select_if(real, is.numeric)
real <- log2(real + 1)

umap_df <- data.frame(UMAP1 = NA, UMAP2 = NA, ID = NA, 
                      groups = NA, datatype = NA, 
                      pilot = NA, model = NA)
for (pilot in selected_pilot) {
  for(model in selected_model){
    cat(paste(pilot, model))
    file <- read.csv(paste("../GeneratedData/", dataname, "_", model, "_", pilot, "_Draw1.csv" , sep = ""), header = F)
    groups <- ifelse(file[,ncol(file)] == 0, g1, g2)
    generated <- file[c(which(groups == g1)[1:n1], which(groups == g2)[1:n2]), 1:ncol(real)]
    colnames(generated) <- colnames(real)
    dat_combine <- rbind(real, generated)
    df <- UMAP_eval(dat_combine, groups = rep(c(g1, g2), c(n1, n2)), log = TRUE, failure = "remove") 
    df <- df$plot_df
    df$pilot <- pilot
    df$model <- model
    umap_df <- as.data.frame(rbind(umap_df, df))
  }
}

umap_df$pilot <- factor(umap_df$pilot, levels = selected_pilot)
umap_df$model <- factor(umap_df$model, levels = selected_model)
umap_df$Data <- paste(umap_df$groups, umap_df$datatype)

umap_df <- as.data.frame(umap_df[-1,])
umap_df$pilotsize <- factor(paste("Pilot size : ", umap_df$pilot, sep = ""), levels = paste("Pilot size : ", selected_pilot, sep = ""))
umapp <- umap_df %>% ggplot(aes(x = UMAP1, y = UMAP2, color = Data))+
  geom_point(size=1)+
  theme(legend.position="bottom")+
  facet_wrap(vars(pilotsize))+
  labs(title="", x = "", y = "")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(legend.position = "bottom")+
  scale_color_brewer(palette = "PiYG")

ggsave(paste(default, dataname, "_UMAP.pdf", sep = ""),
       gg_default(umapp), width = 10, height=5, unit="in")

#### Optimal evaluation vis ####
load("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Interns/2021 Summer/Generative Models/DGMmicroRNA/FinalResAug2023/LAMLPositive_2_Optimal.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Interns/2021 Summer/Generative Models/DGMmicroRNA/FinalResAug2023/PRADPositive_2_Optimal.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Interns/2021 Summer/Generative Models/DGMmicroRNA/FinalResAug2023/BRCAPRADPositive_3_Optimal.RData")
dataname <- "LAMLPositive_2_PRAD_Positive_3"
res_laml <- df_proc(laml_pos_opt)
res_prad <- df_proc(prad_pos_opt)
res <- rbind(res_laml, res_prad)
res$Cancer <- rep(c("LAML", "PRAD"), c(nrow(res_laml), nrow(res_prad)))

dataname <- "BRCAPRADPositive_3"
# change to cancer two types df_proc
res <- df_proc(brcaprad_pos_opt)
res$Cancer <- "BRCAPRAD"

mycolors <- c(brewer.pal(n = 8, name = "Blues")[c(6)], brewer.pal(n = 8, name = "Oranges")[c(6)], brewer.pal(n = 8, name = "Greens")[c(6)])
p <- ggplot(data = res) + 
  geom_violin(aes(x = pilot_size, y = Rate, fill = model), position = position_dodge(0.6))+
  # geom_boxplot(aes(x = pilot_size, y = Rate, fill = model), show.legend = T, width = 0.6)+
  facet_grid(row = vars(Metric), col = vars(Cancer))+
  theme_bw()+
  labs(x = "Pilot size", y = "", switch = "y", fill = "Models")+
  scale_fill_manual(values = mycolors)+
  theme(panel.border = element_rect(colour = "black", fill = NA))+
  theme(legend.position = "bottom")+
  ylim(0,1)
pdf(file = paste("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Interns/2021 Summer/Generative Models/DGMmicroRNA/FinalResAug2023/"
                 , dataname, "_Optimal.pdf", sep = ""), width = 4, height = 9)
print(gg_default(p))
dev.off()



#### Transfer learning vis ####
load("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Interns/2021 Summer/Generative Models/DGMmicroRNA/FinalResAug2023/SKCM_Transfer.RData")
load("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Interns/2021 Summer/Generative Models/DGMmicroRNA/FinalResAug2023/BRCA_Transfer.RData")
skcm_transfer <- skcm_transfer[skcm_transfer$modelname == "VAE1-10",]
brca_transfer <- brca_transfer[brca_transfer$modelname == "VAE1-10",]
res1 <- df_proc(skcm_transfer)
res2 <- df_proc(brca_transfer)
res1$cancer = "miRNA-SKCM"
res2$cancer = "miRNA-BRCA"
res <- rbind(res1, res2)
res$cancer <- factor(res$cancer, levels = c("miRNA-SKCM","miRNA-BRCA"))
dataname <- "miRNA_Transfer"
res$transfer <- factor(res$transfer, levels = c("Original","transfromLAML","transfromLAMLBRCAPRAD","transfromPRAD","transfromSKCMLAMLPRAD"), 
                       labels = c("None", "fromLAML","fromLAMLBRCAPRAD","fromPRAD","fromSKCMLAMLPRAD"))
p <- ggplot(data = res) + 
  geom_boxplot(aes(x = pilot_size, y = Rate, fill = transfer), show.legend = T, width = 0.6)+
  facet_grid(col = vars(Metric), row = vars(cancer))+
  theme_bw()+
  labs(x = "Pilot size", y = "", switch = "y", fill = "Transfer")+
  scale_color_npg()+
  theme(panel.border = element_rect(colour = "black", fill = NA))+
  theme(legend.position = "bottom")
pdf(file = paste("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Interns/2021 Summer/Generative Models/DGMmicroRNA/FinalResAug2023/"
                 , dataname, ".pdf", sep = ""), width = 8, height = 5)
print(gg_default(p))
dev.off()


#RNA transfer
load("/Users/yunhuiqi/Library/Mobile Documents/com~apple~CloudDocs/PhD/Interns/2021 Summer/Generative Models/DGMmicroRNA/Results2023/RNABRCASPLIT_Transfer.RData")
load("/Users/yunhuiqi/Library/Mobile Documents/com~apple~CloudDocs/PhD/Interns/2021 Summer/Generative Models/DGMmicroRNA/Results2023/RNAPRADSPLIT_Transfer.RData")

rnares1 <- df_proc(rnabrca_transfer)
rnares2 <- df_proc(rnaprad_transfer)
rnares1$cancer = "RNABRCA"
rnares2$cancer = "RNAPRAD"
rnares1$transfer <- ifelse(rnares1$transfer== "Original", "None", "fromRNAPRAD")
rnares2$transfer <- ifelse(rnares2$transfer== "Original", "None", "fromRNABRCA")
dataname = "RNA_Transfer"
rnares = rbind(rnares1, rnares2)
rnares$transfer <- factor(rnares$transfer, levels = c("None", "fromRNAPRAD","fromRNABRCA"))
prna <- ggplot(data = rnares) + 
  geom_boxplot(aes(x = pilot_size, y = Rate, fill = transfer), show.legend = T, width = 0.6)+
  facet_grid(col = vars(Metric), row = vars(cancer))+
  theme_bw()+
  labs(x = "Pilot size", y = "", switch = "y", fill = "Transfer")+
  scale_color_npg()+
  theme(panel.border = element_rect(colour = "black", fill = NA))+
  theme(legend.position = "bottom")
pdf(file = paste("~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Interns/2021 Summer/Generative Models/DGMmicroRNA/FinalResAug2023/"
                 , dataname, ".pdf", sep = ""), width = 5, height = 5)
print(gg_default(prna))
dev.off()
#### RNA evaluation vis ####
vis_all_RNA <- function(dataname, nickname, cancer_type = c(1, 2)){
  # This function visualization the results of DGM evaluation
  # dataname = "RNABRCAPositive_5-2", nickname = "rnabrca_pos", cancer_type = 1
  
  # Step1: load data ####
  default <- "~/Library/Mobile Documents/com~apple~CloudDocs/PhD/Interns/2021 Summer/Generative Models/DGMmicroRNA/FinalResAug2023/"
  load(paste(default, dataname, ".RData", sep = ""))
  load(paste(default, dataname, "_Gaussianhead.RData", sep = ""))
  load(paste(default, dataname, "_Normalization.RData", sep = ""))
  overall <- get(nickname)
  offline <- get(paste(nickname, "_Gaussianhead", sep = ""))
  norm <- get(paste(nickname, "_norm", sep = ""))
  
  # Step2: generate dataframes for plotting  ####
  if(cancer_type == 1){
    # Evaluation criteria: 1-Prop.successful.genes, 2-CCC.pos.control, 3-cARI, 4-MAD
    mycolors <- c(brewer.pal(n = 8, name = "Blues")[c(6)], brewer.pal(n = 8, name = "Oranges")[c(6)], brewer.pal(n = 8, name = "Greens")[c(6)])
    plot_width <- c(6, 9, 6, 4)
    plot_height <- c(4, 4, 4, 4)
    df_proc <- function(res){
      res$succ_prop <- 1 - res$fail_prop
      res$pilot_size <- factor(res$pilot_size, levels = as.character(sort(as.numeric(as.character(unique(res$pilot_size))))))
      if("modelname" %in% colnames(res)){model <- res[,"modelname"]}
      else{model <- res[,"model"]}
      lvs <- c("VAE1-100","MAF","WGANGP")
      res$model <- factor(model, levels = lvs)
      res$modeltype <- rep(NA, nrow(res))
      for (i in 1:nrow(res)) {
        if(substr(res$model[i], 1,3)=="VAE"){res$modeltype[i] <- "VAEs"}
        else if(substr(res$model[i], 1,3) %in% c("GAN", "WGA")){res$modeltype[i] <- "GANs"}
        else{res$modeltype[i] <- "FLOWs"}
      }
      res$modeltype <- factor(res$modeltype, levels = c("VAEs","FLOWs","GANs"))
      res$ARI[res$ARI<0 | is.na(res$ARI)] <- 0
      res <- tidyr::pivot_longer(res, cols = c("ARI",  "succ_prop"), names_to = "Metric", values_to = "Rate")
      res$Metric <- factor(res$Metric, levels = c("succ_prop","ARI"),
                           labels = c("1 - Pct(0-genes)", "cARI"))
      if("normalization" %in% colnames(res)){res$normalization <- factor(res$normalization, levels=c('None','TC','TMM','UQ'))}
      return(res)
    }
    heatmap_proc <- function(res){
      if("modelname" %in% colnames(res)){model <- res[,"modelname"]}
      else{model <- res[,"model"]}
      lvs <- c("VAE1-100","MAF","WGANGP")
      res$model <- factor(model, levels = lvs)
      res$pilot_size <- factor(res$pilot_size, levels = as.character(sort(as.numeric(as.character(unique(res$pilot_size))))))
      res$modeltype <- rep(NA, nrow(res))
      for (i in 1:nrow(res)) {
        if(substr(res$model[i], 1,3)=="VAE"){res$modeltype[i] <- "VAEs"}
        else if(substr(res$model[i], 1,3) %in% c("GAN", "WGA")){res$modeltype[i] <- "GANs"}
        else{res$modeltype[i] <- "FLOWs"}
      }
      res$modeltype <- factor(res$modeltype, levels = c("VAEs","FLOWs","GANs"))
      res <- res %>% group_by(pilot_size, model, modeltype) %>%
        summarise(
          Mean = mean(ks_mean,na.rm=TRUE),
          Sd = mean(ks_sd,na.rm=TRUE),
          Sparsity = mean(ks_zero,na.rm=TRUE)
        )
      
      res <- tidyr::pivot_longer(res, cols = c("Mean", "Sd", "Sparsity"), names_to = "Metric", values_to = "Rate")
      res$Metric <- factor(res$Metric, levels = (c("Sparsity","Sd", "Mean")), labels = c("MAD(Sparsity)","MAD(SD)","MAD(Mean)"))
      res$text <- round(res$Rate, 3)
      res$text[res$Rate > 6 | is.na(res$Rate)] <- ">6"
      res$Rate[res$Rate > 6 | is.na(res$Rate)] <- 6
      res$modelorder <- factor(res$model, levels=rev(levels(res$model)))
      res$modeltype <- "All models"
      if("normalization" %in% colnames(res)){res$normalization <- factor(res$normalization, levels=c('None','TC','TMM','UQ','DESeq'))}
      return(res)
    }
  }
  else{
    # Evaluation criteria: 1-Prop.successful.genes, 2-CCC.pos.control, 3-ARI, 4-MAD, 5-ccc.log10pvalue, 6-ccc.log2FC
    plot_width <- c(4,6,4,4)
    plot_height <- c(6,4,6,6)
    mycolors <- c(brewer.pal(n = 8, name = "Blues")[c(6)],brewer.pal(n = 8, name = "Oranges")[6])
    df_proc <- function(res){
      res$succ_prop <- 1 - res$fail_prop
      res$pilot_size <- factor(res$pilot_size, levels = as.character(sort(as.numeric(as.character(unique(res$pilot_size))))))
      if("modelname" %in% colnames(res)){model <- res[,"modelname"]}
      else{model <- res[,"model"]}
      lvs <- c("CVAE1-100","MAF")
      res$model <- factor(model, levels = lvs)
      res$modeltype <- rep(NA, nrow(res))
      for (i in 1:nrow(res)) {
        if(substr(res$model[i], 1,3)=="CVA"){res$modeltype[i] <- "VAEs"}
        else{res$modeltype[i] <- "FLOWs"}
      }
      res$modeltype <- factor(res$modeltype, levels = c("VAEs","FLOWs"))
      if("normalization" %in% colnames(res)){res$normalization <- factor(res$normalization, levels=c('None','TC','UQ','TMM'))}
      res$ARI[res$ARI<0 | is.na(res$ARI)] <- 0
      res$ccc_log10pvalue[res$ccc_log10pvalue<0 | is.na(res$ccc_log10pvalue)] <- 0
      res$ccc_log2FC[res$ccc_log2FC<0 | is.na(res$ccc_log2FC)] <- 0
      res <- tidyr::pivot_longer(res, cols = c("ccc_log10pvalue", "ccc_log2FC","ARI", "succ_prop"), names_to = "Metric", values_to = "Rate")
      res$Metric <- factor(res$Metric, levels = c("succ_prop","ARI","ccc_log10pvalue", "ccc_log2FC"),
                           labels = c("1 - Pct(0-genes)",  "ARI","CCC of -log10(pvalue)", "CCC of log2FC"))
      return(res)
    }
    
    heatmap_proc <- function(res){
      if("modelname" %in% colnames(res)){model <- res[,"modelname"]}
      else{model <- res[,"model"]}
      lvs <- c("CVAE1-100","MAF")
      res$model <- factor(model, levels = lvs)
      res$pilot_size <- factor(res$pilot_size, levels = as.character(sort(as.numeric(as.character(unique(res$pilot_size))))))
      res$modeltype <- rep(NA, nrow(res))
      for (i in 1:nrow(res)) {
        if(substr(res$model[i], 1,3)=="CVA"){res$modeltype[i] <- "VAEs"}
        else{res$modeltype[i] <- "FLOWs"}
      }
      res$modeltype <- factor(res$modeltype, levels = c("VAEs","FLOWs"))
      res <- res %>% group_by(pilot_size, model, modeltype) %>%
        summarise(
          Mean = mean(ks_mean,na.rm=TRUE),
          Sd = mean(ks_sd,na.rm=TRUE),
          Sparsity = mean(ks_zero,na.rm=TRUE)
        )
      
      res <- tidyr::pivot_longer(res, cols = c("Mean", "Sd", "Sparsity"), names_to = "Metric", values_to = "Rate")
      res$Metric <- factor(res$Metric, levels = (c("Sparsity","Sd", "Mean")), labels = c("MAD(Sparsity)","MAD(SD)","MAD(Mean)"))
      res$text <- round(res$Rate, 3)
      res$text[res$Rate > 6 | is.na(res$Rate)] <- ">6"
      res$Rate[res$Rate > 6 | is.na(res$Rate)] <- 6
      if("normalization" %in% colnames(res)){res$normalization <- factor(res$normalization, levels=c('None','TC','UQ','TMM',"DESeq"))}
      return(res)
    }
  }
  # Step 3: make graphs  ####
  gg_default <- function(p){
    pp <- p +
      theme(legend.position = "bottom")+
      theme(strip.text.x = element_text(face = "bold", size = 12),
            strip.text.y = element_text(face = "bold", size = 10),
            legend.title=element_text(size = 10), 
            legend.text=element_text(size = 9),
            axis.text = element_text(size = 12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank())
    return(pp)
  }
  #### p11 overall boxplot ####
  res <- df_proc(overall)
  p11 <- ggplot(data = res[res$pilot_size %in% c("50","150","250"),])+ 
    geom_violin(aes(x = pilot_size, y = Rate, fill = model), show.legend = T, position = position_dodge(0.9))+
    facet_grid(row = vars(Metric), col = vars(modeltype))+
    theme_bw()+
    labs(x = "Pilot size", y = "", switch = "y", fill = "Models")+
    scale_fill_manual(values = mycolors)+
    theme(panel.border = element_rect(colour = "black", fill = NA))+
    guides(fill=guide_legend(nrow = 1,byrow = TRUE))+
    ylim(0,1)
  pdf(file = paste(default, dataname, ".pdf", sep = ""), width = plot_width[1], height = plot_height[1])
  print(gg_default(p11))
  dev.off()
  #### p12 overall heatmap ####
  res <- heatmap_proc(overall)
  p12 <- ggplot(res[res$pilot_size %in% c("50","150","250"),], aes(x = pilot_size , y = Metric, fill = Rate)) + 
    geom_tile(show.legend = T)+
    labs(y = "",x = "Pilot size", fill = "MAD")+
    facet_wrap(vars(model), dir = "h")+
    geom_text(aes(pilot_size, Metric, label = text), color = "black", size = 4) +
    theme_minimal()+
    scale_fill_gradientn(colours = terrain.colors(10))+
    theme(axis.text.y = element_text(angle = 0, hjust = 1, colour = "black", size = 10, face = "bold"))+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank())
  pdf(file = paste(default, dataname, "_MAD.pdf", sep = ""), width = plot_width[2], height = plot_height[2])
  print(gg_default(p12))
  dev.off()
  
  #### p2 offline violine plot####
  res <- df_proc(offline)
  p2 <- ggplot(data = res, aes(x = pilot_size, y = Rate))+
    geom_violin(aes(color = offline), show.legend = T, position = position_dodge(0.9)) +
    geom_boxplot(aes(color = offline), width = 0.15, position = position_dodge(0.9)) +
    facet_grid(cols=vars(model), rows= vars(Metric), scale = "free_y")+
    theme_bw()+
    labs(x="Pilot size", y="", color = "Offline")+
    scale_y_continuous(position = "left")+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.line = element_blank())+
    scale_color_npg()+
    ylim(0,1)
  pdf(file = paste(default, dataname, "_Offline.pdf", sep = ""), width = plot_width[3], height = plot_height[3])
  print(gg_default(p2))
  dev.off()
  
  #### p3 normalization violin plot ####
  res <- df_proc(norm)
  p3 <- ggplot(data = res[res$normalization!="DESeq",])+
    geom_violin(aes(x = pilot_size, y = Rate, fill = normalization ), position = position_dodge(0.9))+
    facet_grid(rows = vars(Metric), cols = vars(model),scales = "free_y")+
    theme_bw()+
    labs(x="Pilot size", y="", fill="Normalization")+
    scale_y_continuous(position = "left")+
    theme(panel.border = element_rect(colour = "black", fill = NA),
          axis.line = element_blank())+
    scale_fill_npg()+
    ylim(0,1)
  pdf(file = paste(default, dataname, "_Normalization.pdf", sep = ""), width = plot_width[4], height = plot_height[4])
  print(gg_default(p3))
  dev.off()
  
}

vis_all_RNA("RNABRCAPositive_5-2", "rnabrca_pos", 1)
vis_all_RNA("RNAPRADPositive_5-2", "rnaprad_pos", 1)
vis_all_RNA("RNABRCAPRADPositive_5-2", "rnabrcaprad_pos", 2)


########################################################################################################################################
# 16. Stage-based FD Z-scores
## This R code calculates Z-scores for stage-level analyses
## it produces supplementary table S7
#######################################################################################################################################

## Import packages
library(tidyverse)
library(writexl)

# Load empirical and null results
load(file = "~/Stage_variation_metrics.RData")
load(file = "~/Stage_variation_long_metrics.RData")
load(file = "~/Mean_Stage_metrics.RData")
load(file = "~/Stage_null_variation_metrics.RData")
load(file = "~/Stage_null_variation_long_metrics.RData")
load(file = "~/Mean_Stage_Null_metrics.RData")

# Calculate Z scores between empirical analyses and null model per Stage
# Form empty dataframe
stage.emp.null<- as.data.frame(matrix(data= NA,nrow= 22, ncol= 16, dimnames= list(c("Danian","Selandian","Thanetian",
                                                                                    "Ypresian","Lutetian","Bartonian","Priabonian",
                                                                                    "Rupelian","Chattian","Aquitanian","Burdigalian",
                                                                                    "Langhian","Serravallian","Tortonian","Messinian",
                                                                                    "Zanclean","Piacenzian","Gelasian","Calabrian",
                                                                                    "Chibanian","Late/Upper","Recent"),
                                                                           c("nb_sp","nb_fe", "fred", "fored", "fvuln", "fric", "fori","fspe",
                                                                             "nb_sp_Z","nb_fe_Z", "fred_Z", "fored_Z", "fvuln_Z", "fric_Z", "fori_Z","fspe_Z"))))
# Separate by Stage for each analysis
## Empirical 
FDmetrics_taxonvar_1<- Res_FDmetrics_StageVar[grep("Danian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_2<- Res_FDmetrics_StageVar[grep("Selandian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_3<- Res_FDmetrics_StageVar[grep("Thanetian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_4<- Res_FDmetrics_StageVar[grep("Ypresian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_5<- Res_FDmetrics_StageVar[grep("Lutetian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_6<- Res_FDmetrics_StageVar[grep("Bartonian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_7<- Res_FDmetrics_StageVar[grep("Priabonian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_8<- Res_FDmetrics_StageVar[grep("Rupelian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_9<- Res_FDmetrics_StageVar[grep("Chattian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_10<- Res_FDmetrics_StageVar[grep("Aquitanian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_11<- Res_FDmetrics_StageVar[grep("Burdigalian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_12<- Res_FDmetrics_StageVar[grep("Langhian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_13<- Res_FDmetrics_StageVar[grep("Serravallian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_14<- Res_FDmetrics_StageVar[grep("Tortonian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_15<- Res_FDmetrics_StageVar[grep("Messinian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_16<- Res_FDmetrics_StageVar[grep("Zanclean", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_17<- Res_FDmetrics_StageVar[grep("Piacenzian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_18<- Res_FDmetrics_StageVar[grep("Gelasian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_19<- Res_FDmetrics_StageVar[grep("Calabrian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_20<- Res_FDmetrics_StageVar[grep("Chibanian", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_21<- Res_FDmetrics_StageVar[grep("Late/Upper", Res_FDmetrics_StageVar$Stage), (2:9)]
FDmetrics_taxonvar_Rec<- Res_FDmetrics_StageVar[grep("Recent", Res_FDmetrics_StageVar$Stage), (2:9)]
## Null model
Res_FDmetrics_StageNull_1<- Res_FDmetrics_StageNull[grep("Danian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_2<- Res_FDmetrics_StageNull[grep("Selandian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_3<- Res_FDmetrics_StageNull[grep("Thanetian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_4<- Res_FDmetrics_StageNull[grep("Ypresian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_5<- Res_FDmetrics_StageNull[grep("Lutetian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_6<- Res_FDmetrics_StageNull[grep("Bartonian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_7<- Res_FDmetrics_StageNull[grep("Priabonian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_8<- Res_FDmetrics_StageNull[grep("Rupelian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_9<- Res_FDmetrics_StageNull[grep("Chattian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_10<- Res_FDmetrics_StageNull[grep("Aquitanian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_11<- Res_FDmetrics_StageNull[grep("Burdigalian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_12<- Res_FDmetrics_StageNull[grep("Langhian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_13<- Res_FDmetrics_StageNull[grep("Serravallian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_14<- Res_FDmetrics_StageNull[grep("Tortonian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_15<- Res_FDmetrics_StageNull[grep("Messinian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_16<- Res_FDmetrics_StageNull[grep("Zanclean", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_17<- Res_FDmetrics_StageNull[grep("Piacenzian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_18<- Res_FDmetrics_StageNull[grep("Gelasian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_19<- Res_FDmetrics_StageNull[grep("Calabrian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_20<- Res_FDmetrics_StageNull[grep("Chibanian", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_21<- Res_FDmetrics_StageNull[grep("Late/Upper", Res_FDmetrics_StageNull$Stage), (2:9)]
Res_FDmetrics_StageNull_Rec<- Res_FDmetrics_StageNull[grep("Recent", Res_FDmetrics_StageNull$Stage), (2:9)]

# Make empirical dataframes - means or medians of empirical & null models
FDemp <- Res_FDmetrics_StageVar %>% 
  group_by(Stage) %>%
  summarise(Sp_med = median(nb_sp),
            FE_med = median(nb_fe),
            Red_med = median(fred),
            Ored_med = median(fored),
            Vul_med = median(fvuln),
            FRic_med = median(fric),
            Fori_med = median(fori),
            Fspe_med = median(fspe))
FDemp <- FDemp %>% 
  column_to_rownames(var = "Stage")

FDnull <- Res_FDmetrics_StageNull %>% 
  group_by(Stage) %>%
  summarise(Sp_med = median(nb_sp),
            FE_med = median(nb_fe),
            Red_med = median(fred),
            Ored_med = median(fored),
            Vul_med = median(fvuln),
            FRic_med = median(fric),
            Fori_med = median(fori),
            Fspe_med = median(fspe))
FDnull <- FDnull %>% 
  column_to_rownames(var = "Stage")

#Calculate empirical differences; use FDemp as Stages are rownames rather than a column
for(e in 1:8) { 
  for (g in 1:22){
    stage.emp.null[g,e]<- as.numeric(FDemp[g,e])-as.numeric(FDnull[g,e])
    
  }
}


#write function to calculate slopes for null model values
slopes_fun<- function(data1, data2){
  output<- matrix(data= NA, nrow= dim(data1)[1], ncol= dim(data1)[2])
  for (k in 1:8) {
    output[,k]<- data2[,k]- data1[,k]
  }
  return(output)
}

N_1_null_Slopes<- slopes_fun(FDmetrics_taxonvar_1, Res_FDmetrics_StageNull_1)
N_2_null_Slopes<- slopes_fun(FDmetrics_taxonvar_2, Res_FDmetrics_StageNull_2)
N_3_null_Slopes<- slopes_fun(FDmetrics_taxonvar_3, Res_FDmetrics_StageNull_3)
N_4_null_Slopes<- slopes_fun(FDmetrics_taxonvar_4, Res_FDmetrics_StageNull_4)
N_5_null_Slopes<- slopes_fun(FDmetrics_taxonvar_5, Res_FDmetrics_StageNull_5)
N_6_null_Slopes<- slopes_fun(FDmetrics_taxonvar_6, Res_FDmetrics_StageNull_6)
N_7_null_Slopes<- slopes_fun(FDmetrics_taxonvar_7, Res_FDmetrics_StageNull_7)
N_8_null_Slopes<- slopes_fun(FDmetrics_taxonvar_8, Res_FDmetrics_StageNull_8)
N_9_null_Slopes<- slopes_fun(FDmetrics_taxonvar_9, Res_FDmetrics_StageNull_9)
N_10_null_Slopes<- slopes_fun(FDmetrics_taxonvar_10, Res_FDmetrics_StageNull_10)
N_11_null_Slopes<- slopes_fun(FDmetrics_taxonvar_11, Res_FDmetrics_StageNull_11)
N_12_null_Slopes<- slopes_fun(FDmetrics_taxonvar_12, Res_FDmetrics_StageNull_12)
N_13_null_Slopes<- slopes_fun(FDmetrics_taxonvar_13, Res_FDmetrics_StageNull_13)
N_14_null_Slopes<- slopes_fun(FDmetrics_taxonvar_14, Res_FDmetrics_StageNull_14)
N_15_null_Slopes<- slopes_fun(FDmetrics_taxonvar_15, Res_FDmetrics_StageNull_15)
N_16_null_Slopes<- slopes_fun(FDmetrics_taxonvar_16, Res_FDmetrics_StageNull_16)
N_17_null_Slopes<- slopes_fun(FDmetrics_taxonvar_17, Res_FDmetrics_StageNull_17)
N_18_null_Slopes<- slopes_fun(FDmetrics_taxonvar_18, Res_FDmetrics_StageNull_18)
N_19_null_Slopes<- slopes_fun(FDmetrics_taxonvar_19, Res_FDmetrics_StageNull_19)
N_20_null_Slopes<- slopes_fun(FDmetrics_taxonvar_20, Res_FDmetrics_StageNull_20)
N_21_null_Slopes<- slopes_fun(FDmetrics_taxonvar_21, Res_FDmetrics_StageNull_21)
N_Rec_null_Slopes<- slopes_fun(FDmetrics_taxonvar_Rec, Res_FDmetrics_StageNull_Rec)

## Calculate Z scores
for (i in 1:8){
  stage.emp.null[1,i+8]<- (stage.emp.null[1,i]- median(N_1_null_Slopes[,i], na.rm = TRUE))/sd(N_1_null_Slopes[,i], na.rm= TRUE)}

for (i in 1:8){
  stage.emp.null[2,i+8]<- (stage.emp.null[2,i]- median(N_2_null_Slopes[,i],  na.rm = TRUE))/sd(N_2_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[3,i+8]<- (stage.emp.null[3,i]- median(N_3_null_Slopes[,i], na.rm = TRUE))/sd(N_3_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[4,i+8]<- (stage.emp.null[4,i]- median(N_4_null_Slopes[,i],  na.rm = TRUE))/sd(N_4_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[5,i+8]<- (stage.emp.null[5,i]- median(N_5_null_Slopes[,i],  na.rm = TRUE))/sd(N_5_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[6,i+8]<- (stage.emp.null[6,i]- median(N_6_null_Slopes[,i],  na.rm = TRUE))/sd(N_6_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[7,i+8]<- (stage.emp.null[7,i]- median(N_7_null_Slopes[,i], na.rm = TRUE))/sd(N_7_null_Slopes[,i], na.rm= TRUE)}

for (i in 1:8){
  stage.emp.null[8,i+8]<- (stage.emp.null[8,i]- median(N_8_null_Slopes[,i],  na.rm = TRUE))/sd(N_8_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[9,i+8]<- (stage.emp.null[9,i]- median(N_9_null_Slopes[,i], na.rm = TRUE))/sd(N_9_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[10,i+8]<- (stage.emp.null[10,i]- median(N_10_null_Slopes[,i],  na.rm = TRUE))/sd(N_10_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[11,i+8]<- (stage.emp.null[11,i]- median(N_11_null_Slopes[,i],  na.rm = TRUE))/sd(N_11_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[12,i+8]<- (stage.emp.null[12,i]- median(N_12_null_Slopes[,i],  na.rm = TRUE))/sd(N_12_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[13,i+8]<- (stage.emp.null[13,i]- median(N_13_null_Slopes[,i], na.rm = TRUE))/sd(N_13_null_Slopes[,i], na.rm= TRUE)}

for (i in 1:8){
  stage.emp.null[14,i+8]<- (stage.emp.null[14,i]- median(N_14_null_Slopes[,i],  na.rm = TRUE))/sd(N_14_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[15,i+8]<- (stage.emp.null[15,i]- median(N_15_null_Slopes[,i], na.rm = TRUE))/sd(N_15_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[16,i+8]<- (stage.emp.null[16,i]- median(N_16_null_Slopes[,i],  na.rm = TRUE))/sd(N_16_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[17,i+8]<- (stage.emp.null[17,i]- median(N_17_null_Slopes[,i],  na.rm = TRUE))/sd(N_17_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[18,i+8]<- (stage.emp.null[18,i]- median(N_18_null_Slopes[,i],  na.rm = TRUE))/sd(N_18_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[19,i+8]<- (stage.emp.null[19,i]- median(N_19_null_Slopes[,i],  na.rm = TRUE))/sd(N_19_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[20,i+8]<- (stage.emp.null[20,i]- median(N_20_null_Slopes[,i],  na.rm = TRUE))/sd(N_20_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[21,i+8]<- (stage.emp.null[21,i]- median(N_21_null_Slopes[,i],  na.rm = TRUE))/sd(N_21_null_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  stage.emp.null[22,i+8]<- (stage.emp.null[22,i]- median(N_Rec_null_Slopes[,i],  na.rm = TRUE))/sd(N_Rec_null_Slopes[,i],  na.rm = TRUE)}

stage.emp.null

# Extract results as excel worksheet - produces Tables 2 and S4
save(stage.emp.null, file = "~/Stage_FD Z scores.Rdata")
write_xlsx(stage.emp.null, "~/Stage_FD Z scores.xlsx")

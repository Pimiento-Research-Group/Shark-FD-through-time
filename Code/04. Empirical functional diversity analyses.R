########################################################################################################################################
# 04. Empirical functional diversity analyses
## This R code provides functional diversity analyses on our recorded data
## it produces the Rdata files needed to load empirical results for comparison to other models in Figure 2
#######################################################################################################################################

## Import packages
library(tidyverse)
library(tibble)
library(scales)
library(cowplot)
library(doBy)
library(reshape2)
library(RColorBrewer)
library(mFD)
library(ggsci)
library(data.table)
library(foreach)
library(parallel)
library(doParallel)
library(writexl)

# Load data
load(file="~/Cleaned data.RData")

## Form dataframe of traits and occurrences; to be split into two input matrices in the loop
# Occurrence matrix
# Form FTU unique combinations, with epoch included
abun <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,CE,LC,XO,LO,Epoch,Epoch_earliest,Epoch_latest,Current_status)

# Order taxa alphabetically and order Epochs
abun <- abun[order(abun$Taxon_corrected),]
abun$Epoch <- ordered(abun$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                           "Pliocene","Pleistocene", "Recent","nonmatch"))

# Form FTUs column
abun$FTU<-paste(abun$Taxon_corrected, abun$CH, abun$CW, abun$CE, abun$LC, abun$XO, abun$LO, sep="+")

# Cast into wide-format data
wp.abun <- abun %>% 
  pivot_wider(names_from = Epoch, values_from = FTU)

wp.abun <- wp.abun[, c("Taxon_corrected","Epoch_earliest","Epoch_latest","Current_status",
                       "Palaeocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene","Recent","nonmatch")]


# Mark 0 and 1 for absence and presence data in all epochs
wp.abun$Palaeocene[is.na(wp.abun$Palaeocene)] <- "0" 
wp.abun$Eocene[is.na(wp.abun$Eocene)] <- "0"
wp.abun$Oligocene[is.na(wp.abun$Oligocene)] <- "0"
wp.abun$Miocene[is.na(wp.abun$Miocene)] <- "0"
wp.abun$Pliocene[is.na(wp.abun$Pliocene)] <- "0"
wp.abun$Pleistocene[is.na(wp.abun$Pleistocene)] <- "0"
wp.abun$Recent[is.na(wp.abun$Recent)] <- "0"
wp.abun$nonmatch[is.na(wp.abun$nonmatch)] <- "0"

wp.abun <- wp.abun %>% 
  mutate(Palaeocene = replace(Palaeocene,Palaeocene!="0", "1")) %>%
  mutate(Eocene = replace(Eocene,Eocene!="0", "1")) %>% 
  mutate(Oligocene = replace(Oligocene,Oligocene!="0", "1")) %>% 
  mutate(Miocene = replace(Miocene,Miocene!="0", "1")) %>%  
  mutate(Pliocene = replace(Pliocene,Pliocene!="0", "1")) %>%
  mutate(Pleistocene = replace(Pleistocene,Pleistocene!="0", "1")) %>% 
  mutate(Recent = replace(Recent,Recent!="0", "1")) %>%
  mutate(nonmatch = replace(nonmatch,nonmatch!="0", "1"))

## Fill in non-matches
wp.abun <- wp.abun %>% 
  rowwise() %>% 
  mutate(Palaeocene = case_when(Epoch_earliest=="Palaeocene" && nonmatch=="1"~"1",TRUE~Palaeocene)) %>% 
  mutate(Eocene = case_when(Epoch_latest=="Eocene" && nonmatch=="1"~"1",TRUE~Eocene)) %>% 
  mutate(Eocene = case_when(Epoch_earliest=="Eocene" && nonmatch=="1"~"1",TRUE~Eocene)) %>%
  mutate(Oligocene = case_when(Epoch_latest=="Oligocene" && nonmatch=="1"~"1",TRUE~Oligocene)) %>% 
  mutate(Oligocene = case_when(Epoch_earliest=="Oligocene" && nonmatch=="1"~"1",TRUE~Oligocene)) %>% 
  mutate(Miocene = case_when(Epoch_latest=="Miocene" && nonmatch=="1"~"1",TRUE~Miocene)) %>% 
  mutate(Miocene = case_when(Epoch_earliest=="Miocene" && nonmatch=="1"~"1",TRUE~Miocene)) %>% 
  mutate(Pliocene = case_when(Epoch_latest=="Pliocene" && nonmatch=="1"~"1",TRUE~Pliocene))

# Range through via "extant" status
# Automatically range-through extant taxa
wp.abun <- wp.abun %>% 
  rowwise() %>% 
  mutate(Eocene = case_when(Palaeocene == "1" && Current_status == "Extant"~"1", TRUE ~ Eocene)) %>% 
  mutate(Oligocene = case_when(Eocene == "1" && Current_status == "Extant"~"1", TRUE ~ Oligocene)) %>% 
  mutate(Miocene = case_when(Oligocene == "1" && Current_status == "Extant"~"1", TRUE ~ Miocene)) %>% 
  mutate(Pliocene = case_when(Miocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pliocene)) %>% 
  mutate(Pleistocene = case_when(Pliocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pleistocene)) %>% 
  mutate(Recent = case_when(Pleistocene == "1" && Current_status == "Extant"~"1", TRUE ~ Recent))

baskets.range <- wp.abun %>% 
  select(-c(Epoch_earliest,Epoch_latest,nonmatch))

# Remove all genus-level taxa from Recent where all species are accounted for (i.e., make 0)
baskets.range <- baskets.range %>% 
  rowwise() %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Alopias sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Carcharias sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Carcharodon sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Dalatias sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Galeocerdo sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Galeorhinus sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Galeus sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Hemipristis sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Heptranchias sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Isurus sp."~"0", TRUE ~ Recent)) %>%
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Megachasma sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Nebrius sp."~"0", TRUE ~ Recent)) %>% 
  mutate(Recent = case_when(Recent == "1" && Taxon_corrected == "Notorynchus sp."~"0", TRUE ~ Recent))

# Make 0 NAs, merge duplicates, reform NAs as 0
baskets.range[baskets.range == 0] <- NA
coalesce_all_columns <- function(df) {
  return(coalesce(!!! as.list(df)))
}

baskets.range <- baskets.range %>%
  group_by(Taxon_corrected) %>%
  summarise_all(coalesce_all_columns)

baskets.range$Palaeocene[is.na(baskets.range$Palaeocene)] <- "0" 
baskets.range$Eocene[is.na(baskets.range$Eocene)] <- "0"
baskets.range$Oligocene[is.na(baskets.range$Oligocene)] <- "0"
baskets.range$Miocene[is.na(baskets.range$Miocene)] <- "0"
baskets.range$Pliocene[is.na(baskets.range$Pliocene)] <- "0"
baskets.range$Pleistocene[is.na(baskets.range$Pleistocene)] <- "0"
baskets.range$Recent[is.na(baskets.range$Recent)] <- "0"

# Fill in any gaps in range automatically
baskets.range.taxon <- baskets.range %>% 
  rowwise() %>% 
  mutate(Eocene = case_when(Palaeocene == "1" && Miocene == "1"~ "1", TRUE ~ Eocene)) %>%
  mutate(Oligocene = case_when(Palaeocene == "1" && Miocene == "1"~ "1", TRUE ~ Oligocene)) %>%
  mutate(Oligocene = case_when(Eocene == "1" && Miocene == "1"~ "1", TRUE ~ Oligocene)) %>%
  mutate(Miocene = case_when(Oligocene == "1" && Pliocene == "1"~ "1", TRUE ~ Miocene)) %>%
  mutate(Pliocene = case_when(Miocene == "1" && Pleistocene == "1"~ "1", TRUE ~ Pliocene)) %>%
  mutate(Oligocene = case_when(Eocene == "1" && Pliocene == "1"~ "1", TRUE ~ Oligocene)) %>%
  mutate(Miocene = case_when(Eocene == "1" && Pliocene == "1"~ "1", TRUE ~ Miocene))

# Mark and remove species that only occur in Recent
selected_taxa <- baskets.range.taxon %>%
  filter(
    Recent == "1" &
      Palaeocene == "0" &
      Eocene == "0" &
      Oligocene == "0" &
      Miocene == "0" &
      Pliocene == "0" &
      Pleistocene == "0"
  ) %>%
  select(Taxon_corrected)
# Filter taxa that are NOT present in "selected_taxa"
filtered_baskets.range <- baskets.range.taxon %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

# Make taxa row names & tidy dataframe to form final matrix
baskets <- filtered_baskets.range %>% 
  remove_rownames %>% 
  column_to_rownames(var="Taxon_corrected") %>% 
  select(-Current_status)

# Make species richness object (as this will be constant regardless of FTU variation) to add later
baskets.epoch <- filtered_baskets.range %>% 
  select(c(Taxon_corrected,Palaeocene,Eocene,Oligocene,Miocene,Pliocene,Pleistocene,Recent)) %>% 
  rowwise() %>% 
  mutate(Palaeocene = replace(Palaeocene,Palaeocene=="1",Taxon_corrected)) %>% 
  mutate(Eocene = replace(Eocene,Eocene=="1",Taxon_corrected)) %>% 
  mutate(Oligocene = replace(Oligocene,Oligocene=="1",Taxon_corrected)) %>% 
  mutate(Miocene = replace(Miocene,Miocene=="1",Taxon_corrected)) %>% 
  mutate(Pliocene = replace(Pliocene,Pliocene=="1",Taxon_corrected)) %>% 
  mutate(Pleistocene = replace(Pleistocene,Pleistocene=="1",Taxon_corrected)) %>% 
  mutate(Recent = replace(Recent,Recent=="1",Taxon_corrected))

baskets.epoch[baskets.epoch == 0]<-NA

# Form functional taxonomic units (FTUs), grouping by species to account for intraspecific variability
FTU <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,CE,LC,XO,LO)

# Order FTUs alphabetically to help ensure matches between species-trait matrix and occurrence matrix later
FTU <- FTU[order(FTU$Taxon_corrected),]

# Filter taxa that are NOT present in "selected_taxa"
FTU_filtered <- FTU %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

## Run iterations & calculate mean by taxon
# Set up model to run on n iterations
registerDoParallel(cores = 10)
getDoParWorkers()
res.TaxonVar<-NULL
res.TaxonVar<- lapply(1:1000, function(x){
  # Form sliced species-trait matrix
  sharks_traits.var <- FTU_filtered %>%
    group_by(Taxon_corrected) %>%
    slice_sample(n=1, replace = FALSE) %>% 
    remove_rownames %>% 
    column_to_rownames(var="Taxon_corrected")
  
  # Transpose to form final occurrence matrix
  baskets_sharks_weights.var <- t(baskets)
  
  # Make matrix and numeric
  baskets_sharks_weights.var <- data.matrix(baskets_sharks_weights.var, rownames.force = NA)
  class(baskets_sharks_weights.var) <- "numeric"
  
  # Summarise & extract occurrences
  asb_sp_sharks_summ.var <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights.var)
  
  asb_sp_sharks_occ.var <- asb_sp_sharks_summ.var$asb_sp_occ
  
  ############################################################################################
  # Make trait_cat matrix
  sharks_traits_cat.var <- tibble(trait_name = c("CH", "CW", "CE", "LC", "XO", "LO"),
                                  trait_type = c("O", "O", "N", "N", "N", "N"),
                                  trait_weight = c(0.5,0.5,1,1,0.33,0.67))
  
  ############################################################################################
  # Form the functional entities
  sp_to_fe_sharks.var <- mFD::sp.to.fe(
    sp_tr       = sharks_traits.var, 
    tr_cat      = sharks_traits_cat.var, 
    fe_nm_type  = "fe_rank", 
    check_input = TRUE)
  
  # Calculate and display metrics
  alpha_fd_fe_sharks.var <- alpha.fd.fe(asb_sp_sharks_occ.var,
                                        sp_to_fe_sharks.var,
                                        ind_nm = c("fred", "fored", "fvuln"),
                                        check_input = TRUE,
                                        details_returned = TRUE)
  
  # Form dataframe
  FEmetrics.var <- as.data.frame(alpha_fd_fe_sharks.var$asb_fdfe) %>% 
    tibble::rownames_to_column("Epoch")
  FEmetrics.var$Epoch <- ordered(FEmetrics.var$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                               "Pliocene","Pleistocene","Recent"))
  
  # Construct trait distance matrix
  sp_dist_sharks.var <- mFD::funct.dist(
    sp_tr         = sharks_traits.var,
    tr_cat        = sharks_traits_cat.var,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  # Assess quality
  fspaces_quality_sharks.var <- mFD::quality.fspaces(
    sp_dist             = sp_dist_sharks.var,
    maxdim_pcoa         = 10,
    deviation_weighting = c("absolute","squared"),
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  # Return coordinates of each axis
  sp_faxes_coord_sharks.var <- fspaces_quality_sharks.var$"details_fspaces"$"sp_pc_coord"
  
  # Form alpha diversity indices
  alpha_fd_indices_sharks.var <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_faxes_coord_sharks.var[ , c("PC1", "PC2", "PC3")],
    asb_sp_w         = baskets_sharks_weights.var,
    ind_vect         = c("fric", "fori", "fspe"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  fd_ind_values_sharks.var <- alpha_fd_indices_sharks.var$"functional_diversity_indices"
  
  ## Form dataframe
  FD.var <- as.data.frame(fd_ind_values_sharks.var) %>% 
    tibble::rownames_to_column("Epoch")
  FD.var$Epoch <- ordered(FD.var$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                 "Pliocene","Pleistocene","Recent"))
  
  # Form list to merge datasets
  FDind.var = list(FEmetrics.var,FD.var)
  
  # Output
  FDind.var %>% 
    reduce(inner_join, by = "Epoch")
  
  #close loop
})

stopImplicitCluster()

# Merge lists into 1 dataframe

res_df_TaxonVar <- res.TaxonVar %>% 
  bind_rows()

# Remove sp_richn, which is the same as nb_sp - redundant
res_df_TaxonVar$sp_richn<-NULL

# Format dataframe to be loaded for comparison plots
Res_FDmetrics_TaxonVar<- res_df_TaxonVar %>% 
  select(Epoch:fspe)
save(Res_FDmetrics_TaxonVar, file = "~/Taxon_variation_metrics.RData")

# Melt data
FDmetrics_long_TaxonVar<- melt(Res_FDmetrics_TaxonVar, id.vars= "Epoch")
save(FDmetrics_long_TaxonVar, file = "~/Taxon_variation_long_metrics.RData")

# Form dataframe of mean, median and standard deviation of all FD metrics
## Median and SD values produce Table 1
Taxon_var <- Res_FDmetrics_TaxonVar %>% 
  group_by(Epoch) %>%
  summarise(Sp_mean = mean(nb_sp),
            Sp_med = median(nb_sp),
            Sp_sd = sd(nb_sp),
            FE_mean = mean(nb_fe),
            FE_med = median(nb_fe),
            FE_sd = sd(nb_fe),
            Red_mean = mean(fred),
            Red_med = median(fred),
            Red_sd = sd(fred),
            Ored_mean = mean(fored),
            Ored_med = median(fored),
            Ored_sd = sd(fored),
            Vul_mean = mean(fvuln),
            Vul_med = median(fvuln),
            Vul_sd = sd(fvuln),
            FRic_mean = mean(fric),
            FRic_med = median(fric),
            FRic_sd = sd(fric),
            Fori_mean = mean(fori),
            Fori_med = median(fori),
            Fori_sd = sd(fori),
            Fspe_mean = mean(fspe),
            Fspe_med = median(fspe),
            Fspe_sd = sd(fspe))

Taxon_var$Epoch <- ordered(Taxon_var$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                     "Pliocene","Pleistocene","Recent"))

# Save iteration data
save(Taxon_var, file = "~/Median_Taxon_metrics.RData")

# Isolate individual functional diversity metrics
FDmetrics_FE <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "nb_fe")
FDmetrics_Fred <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fred")
FDmetrics_Fored <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fored")
FDmetrics_Fvuln <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fvuln")
FDmetrics_FRic <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fric")
FDmetrics_Fori <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fori")
FDmetrics_Fspe <- FDmetrics_long_TaxonVar %>% 
  filter(variable == "fspe")

# Test data for normal distribution - all show non-normal distribution (P < 0.05)
## Code below focuses on FE, FRed & FRic but can be applied to all metrics; all metrics are assessed per epoch due to high sample size for shapiro wilk test
FE_Taxonvar_pal <- FDmetrics_FE %>% 
  filter(Epoch == "Palaeocene")
shapiro.test(FE_Taxonvar_pal$value)
FRed_Taxonvar_pal <- FDmetrics_Fred %>% 
  filter(Epoch == "Palaeocene")
shapiro.test(FRed_Taxonvar_pal$value)
FRic_Taxonvar_pal <- FDmetrics_FRic %>% 
  filter(Epoch == "Palaeocene")
shapiro.test(FRic_Taxonvar_pal$value)

FE_Taxonvar_eo <- FDmetrics_FE %>% 
  filter(Epoch == "Eocene")
shapiro.test(FE_Taxonvar_eo$value)
FRed_Taxonvar_eo <- FDmetrics_Fred %>% 
  filter(Epoch == "Eocene")
shapiro.test(FRed_Taxonvar_eo$value)
FRic_Taxonvar_eo <- FDmetrics_FRic %>% 
  filter(Epoch == "Eocene")
shapiro.test(FRic_Taxonvar_eo$value)

FE_Taxonvar_oli <- FDmetrics_FE %>% 
  filter(Epoch == "Oligocene")
shapiro.test(FE_Taxonvar_oli$value)
FRed_Taxonvar_oli <- FDmetrics_Fred %>% 
  filter(Epoch == "Oligocene")
shapiro.test(FRed_Taxonvar_oli$value)
FRic_Taxonvar_oli <- FDmetrics_FRic %>% 
  filter(Epoch == "Oligocene")
shapiro.test(FRic_Taxonvar_oli$value)

FE_Taxonvar_mio <- FDmetrics_FE %>% 
  filter(Epoch == "Miocene")
shapiro.test(FE_Taxonvar_mio$value)
FRed_Taxonvar_mio <- FDmetrics_Fred %>% 
  filter(Epoch == "Miocene")
shapiro.test(FRed_Taxonvar_mio$value)
FRic_Taxonvar_mio <- FDmetrics_FRic %>% 
  filter(Epoch == "Miocene")
shapiro.test(FRic_Taxonvar_mio$value)

FE_Taxonvar_plio <- FDmetrics_FE %>% 
  filter(Epoch == "Pliocene")
shapiro.test(FE_Taxonvar_plio$value)
FRed_Taxonvar_plio <- FDmetrics_Fred %>% 
  filter(Epoch == "Pliocene")
shapiro.test(FRed_Taxonvar_plio$value)
FRic_Taxonvar_plio <- FDmetrics_FRic %>% 
  filter(Epoch == "Pliocene")
shapiro.test(FRic_Taxonvar_plio$value)

FE_Taxonvar_ple <- FDmetrics_FE %>% 
  filter(Epoch == "Pleistocene")
shapiro.test(FE_Taxonvar_ple$value)
FRed_Taxonvar_ple <- FDmetrics_Fred %>% 
  filter(Epoch == "Pleistocene")
shapiro.test(FRed_Taxonvar_ple$value)
FRic_Taxonvar_ple <- FDmetrics_FRic %>% 
  filter(Epoch == "Pleistocene")
shapiro.test(FRic_Taxonvar_ple$value)

FE_Taxonvar_rec <- FDmetrics_FE %>% 
  filter(Epoch == "Recent")
shapiro.test(FE_Taxonvar_rec$value)
FRed_Taxonvar_rec <- FDmetrics_Fred %>% 
  filter(Epoch == "Recent")
shapiro.test(FRed_Taxonvar_rec$value)
FRic_Taxonvar_rec <- FDmetrics_FRic %>% 
  filter(Epoch == "Recent")
shapiro.test(FRic_Taxonvar_rec$value)

## Calculate changes in FD across epochs - produces proportional changes in Table 1
slopes.emp<- as.data.frame(matrix(data= NA,nrow= 6, ncol= 16, dimnames= list(c("Palaeocene-Eocene", "Eocene-Oligocene", "Oligocene-Miocene", "Miocene-Pliocene", "Pliocene-Pleistocene","Pleistocene-Recent"),
                                                                             c("nb_sp","nb_fe", "fred", "fored", "fvuln", "fric", "fori","fspe",
                                                                               "nb_sp_Z","nb_fe_Z", "fred_Z", "fored_Z", "fvuln_Z", "fric_Z", "fori_Z","fspe_Z"))))
# Make empirical dataframes - means or medians of empirical & null models
FDemp <- Res_FDmetrics_TaxonVar %>% 
  group_by(Epoch) %>%
  summarise(Sp_med = median(nb_sp),
            FE_med = median(nb_fe),
            Red_med = median(fred),
            Ored_med = median(fored),
            Vul_med = median(fvuln),
            FRic_med = median(fric),
            Fori_med = median(fori),
            Fspe_med = median(fspe))
FDemp <- FDemp %>% 
  column_to_rownames(var = "Epoch")

# Separate by epoch
FDmetrics_taxonvar_Pal<- Res_FDmetrics_TaxonVar[grep("Palaeocene", Res_FDmetrics_TaxonVar$Epoch), (2:9)]
FDmetrics_taxonvar_Eo<- Res_FDmetrics_TaxonVar[grep("Eocene", Res_FDmetrics_TaxonVar$Epoch), (2:9)]
FDmetrics_taxonvar_Oli<- Res_FDmetrics_TaxonVar[grep("Oligocene", Res_FDmetrics_TaxonVar$Epoch), (2:9)]
FDmetrics_taxonvar_Mio<- Res_FDmetrics_TaxonVar[grep("Miocene", Res_FDmetrics_TaxonVar$Epoch), (2:9)]
FDmetrics_taxonvar_Plio<- Res_FDmetrics_TaxonVar[grep("Pliocene", Res_FDmetrics_TaxonVar$Epoch), (2:9)]
FDmetrics_taxonvar_Ple<- Res_FDmetrics_TaxonVar[grep("Pleistocene", Res_FDmetrics_TaxonVar$Epoch), (2:9)]
FDmetrics_taxonvar_Rec<- Res_FDmetrics_TaxonVar[grep("Recent", Res_FDmetrics_TaxonVar$Epoch), (2:9)]

# write "slopes" function (Hedberg et al. 2021) to calculate changes 
slopes_fun<- function(data1, data2){
  output<- matrix(data= NA, nrow= dim(data1)[1], ncol= dim(data1)[2])
  for (k in 1:8) {
    output[,k]<- data2[,k]- data1[,k]
  }
  return(output)
}

# Calculate empirical differences; use FDemp as epochs are rownames rather than a column
for(e in 1:8) { 
  for (g in 1:6){
    slopes.emp[g,e]<- as.numeric(FDemp[g+1,e])-as.numeric(FDemp[g,e])
    
  }
}

N_Pal_Eo_Slopes<- slopes_fun(FDmetrics_taxonvar_Pal, FDmetrics_taxonvar_Eo)
N_Eo_Oli_Slopes<- slopes_fun(FDmetrics_taxonvar_Eo, FDmetrics_taxonvar_Oli)
N_Oli_Mio_Slopes<- slopes_fun(FDmetrics_taxonvar_Oli, FDmetrics_taxonvar_Mio)
N_Mio_Plio_Slopes<- slopes_fun(FDmetrics_taxonvar_Mio, FDmetrics_taxonvar_Plio)
N_Plio_Ple_Slopes<- slopes_fun(FDmetrics_taxonvar_Plio, FDmetrics_taxonvar_Ple)
N_Ple_Rec_Slopes<- slopes_fun(FDmetrics_taxonvar_Ple, FDmetrics_taxonvar_Rec)

## Calculate Z scores
for (i in 1:8){
  slopes.emp[1,i+8]<- (slopes.emp[1,i]- median(N_Pal_Eo_Slopes[,i], na.rm = TRUE))/sd(N_Pal_Eo_Slopes[,i], na.rm= TRUE)}

for (i in 1:8){
  slopes.emp[2,i+8]<- (slopes.emp[2,i]- median(N_Eo_Oli_Slopes[,i],  na.rm = TRUE))/sd(N_Eo_Oli_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  slopes.emp[3,i+8]<- (slopes.emp[3,i]- median(N_Oli_Mio_Slopes[,i], na.rm = TRUE))/sd(N_Oli_Mio_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  slopes.emp[4,i+8]<- (slopes.emp[4,i]- median(N_Mio_Plio_Slopes[,i],  na.rm = TRUE))/sd(N_Mio_Plio_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  slopes.emp[5,i+8]<- (slopes.emp[5,i]- median(N_Plio_Ple_Slopes[,i],  na.rm = TRUE))/sd(N_Plio_Ple_Slopes[,i],  na.rm = TRUE)}

for (i in 1:8){
  slopes.emp[6,i+8]<- (slopes.emp[6,i]- median(N_Ple_Rec_Slopes[,i],  na.rm = TRUE))/sd(N_Ple_Rec_Slopes[,i],  na.rm = TRUE)}

slopes.emp
write_xlsx(slopes.emp, "~/Proportional_changes.xlsx")

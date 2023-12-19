###################################################################################################################
# 11. Sensitivity tests 
## This R code provides the sensitivity tests of number of traits for the shark functional diversity analyses
## it produces Figure S9
###################################################################################################################

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
library(deeptime)

# Load data
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Cleaned data.RData")

# Load empirical analyses
load(file="C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Mean_Taxon_metrics.RData")

# Analyses without CH
# Occurrence matrix
# Form FTU unique combinations, with epoch included
abun_CH <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CW,CE,LC,XO,LO,Epoch,Epoch_earliest,Epoch_latest,Current_status)

# Order taxa alphabetically like before and order Epochs
abun_CH <- abun_CH[order(abun_CH$Taxon_corrected),]
abun_CH$Epoch <- ordered(abun_CH$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                 "Pliocene","Pleistocene", "Recent","nonmatch"))

# Form FTUs column
abun_CH$FTU<-paste(abun_CH$Taxon_corrected, abun_CH$CW, abun_CH$CE, abun_CH$LC, abun_CH$XO, abun_CH$LO, sep="+")

# Cast into wide-format data
wp.abun_CH <- abun_CH %>% 
  pivot_wider(names_from = Epoch, values_from = FTU)

wp.abun_CH <- wp.abun_CH[, c("Taxon_corrected","Epoch_earliest","Epoch_latest","Current_status",
                             "Palaeocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene","Recent","nonmatch")]


# Mark 0 and 1 for absence and presence data in all epochs
wp.abun_CH$Palaeocene[is.na(wp.abun_CH$Palaeocene)] <- "0" 
wp.abun_CH$Eocene[is.na(wp.abun_CH$Eocene)] <- "0"
wp.abun_CH$Oligocene[is.na(wp.abun_CH$Oligocene)] <- "0"
wp.abun_CH$Miocene[is.na(wp.abun_CH$Miocene)] <- "0"
wp.abun_CH$Pliocene[is.na(wp.abun_CH$Pliocene)] <- "0"
wp.abun_CH$Pleistocene[is.na(wp.abun_CH$Pleistocene)] <- "0"
wp.abun_CH$Recent[is.na(wp.abun_CH$Recent)] <- "0"
wp.abun_CH$nonmatch[is.na(wp.abun_CH$nonmatch)] <- "0"

wp.abun_CH <- wp.abun_CH %>% 
  mutate(Palaeocene = replace(Palaeocene,Palaeocene!="0", "1")) %>%
  mutate(Eocene = replace(Eocene,Eocene!="0", "1")) %>% 
  mutate(Oligocene = replace(Oligocene,Oligocene!="0", "1")) %>% 
  mutate(Miocene = replace(Miocene,Miocene!="0", "1")) %>%  
  mutate(Pliocene = replace(Pliocene,Pliocene!="0", "1")) %>%
  mutate(Pleistocene = replace(Pleistocene,Pleistocene!="0", "1")) %>% 
  mutate(Recent = replace(Recent,Recent!="0", "1")) %>%
  mutate(nonmatch = replace(nonmatch,nonmatch!="0", "1"))

## Fill in non-matches
wp.abun_CH <- wp.abun_CH %>% 
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
wp.abun_CH <- wp.abun_CH %>% 
  rowwise() %>% 
  mutate(Eocene = case_when(Palaeocene == "1" && Current_status == "Extant"~"1", TRUE ~ Eocene)) %>% 
  mutate(Oligocene = case_when(Eocene == "1" && Current_status == "Extant"~"1", TRUE ~ Oligocene)) %>% 
  mutate(Miocene = case_when(Oligocene == "1" && Current_status == "Extant"~"1", TRUE ~ Miocene)) %>% 
  mutate(Pliocene = case_when(Miocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pliocene)) %>% 
  mutate(Pleistocene = case_when(Pliocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pleistocene)) %>% 
  mutate(Recent = case_when(Pleistocene == "1" && Current_status == "Extant"~"1", TRUE ~ Recent))

baskets.range <- wp.abun_CH %>% 
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

Spp_richness <- c(length(unique(baskets.epoch$Palaeocene)),
                  length(unique(baskets.epoch$Eocene)),
                  length(unique(baskets.epoch$Oligocene)),
                  length(unique(baskets.epoch$Miocene)),
                  length(unique(baskets.epoch$Pliocene)),
                  length(unique(baskets.epoch$Pleistocene)),
                  length(unique(baskets.epoch$Recent)))

# Form functional taxonomic units (FTUs), grouping by species to account for intraspecific variability
FTU_CH <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CW,CE,LC,XO,LO)

# Order FTUs alphabetically to help ensure matches between species-trait matrix and occurrence matrix later
FTU_CH <- FTU_CH[order(FTU_CH$Taxon_corrected),]

# Filter taxa that are NOT present in "selected_taxa"
FTU_filtered <- FTU_CH %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

## Run iterations & calculate mean by taxon
# Set up model to run on n iterations
registerDoParallel(cores = 10)
getDoParWorkers()
res.CH<-NULL
res.CH<- lapply(1:1000, function(x){
  # Form sliced species-trait matrix
  sharks_traits.CH <- FTU_filtered %>%
    group_by(Taxon_corrected) %>%
    slice_sample(n=1, replace = FALSE) %>% 
    remove_rownames %>% 
    column_to_rownames(var="Taxon_corrected")
  
  # Transpose to form final occurrence matrix
  baskets_sharks_weights.CH <- t(baskets)
  
  # Make matrix and numeric
  baskets_sharks_weights.CH <- data.matrix(baskets_sharks_weights.CH, rownames.force = NA)
  class(baskets_sharks_weights.CH) <- "numeric"
  
  # Summarise & extract occurrences
  asb_sp_sharks_summ.CH <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights.CH)
  
  asb_sp_sharks_occ.CH <- asb_sp_sharks_summ.CH$asb_sp_occ
  
  ############################################################################################
  # Make trait_cat matrix
  sharks_traits_cat.CH <- tibble(trait_name = c("CW", "CE", "LC", "XO", "LO"),
                                 trait_type = c("O", "N", "N", "N", "N"),
                                 trait_weight = c(0.5,1,1,0.33,0.67))
  
  ############################################################################################
  # Form the functional entities
  sp_to_fe_sharks.CH <- mFD::sp.to.fe(
    sp_tr       = sharks_traits.CH, 
    tr_cat      = sharks_traits_cat.CH, 
    fe_nm_type  = "fe_rank", 
    check_input = TRUE)
  
  # Calculate and display metrics
  alpha_fd_fe_sharks.CH <- alpha.fd.fe(asb_sp_sharks_occ.CH,
                                       sp_to_fe_sharks.CH,
                                       ind_nm = c("fred", "fored", "fvuln"),
                                       check_input = TRUE,
                                       details_returned = TRUE)
  
  # Form dataframe
  FEmetrics.CH <- as.data.frame(alpha_fd_fe_sharks.CH$asb_fdfe) %>% 
    tibble::rownames_to_column("Epoch")
  FEmetrics.CH$Epoch <- ordered(FEmetrics.CH$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                             "Pliocene","Pleistocene","Recent"))
  
  # Construct trait distance matrix
  sp_dist_sharks.CH <- mFD::funct.dist(
    sp_tr         = sharks_traits.CH,
    tr_cat        = sharks_traits_cat.CH,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  # Assess quality
  fspaces_quality_sharks.CH <- mFD::quality.fspaces(
    sp_dist             = sp_dist_sharks.CH,
    maxdim_pcoa         = 10,
    deviation_weighting = c("absolute","squared"),
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  # Return coordinates of each PC/axis
  sp_faxes_coord_sharks.CH <- fspaces_quality_sharks.CH$"details_fspaces"$"sp_pc_coord"
  
  # Form alpha diversity indices
  alpha_fd_indices_sharks.CH <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_faxes_coord_sharks.CH[ , c("PC1", "PC2", "PC3")],
    asb_sp_w         = baskets_sharks_weights.CH,
    ind_vect         = c("fric", "fori", "fspe"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  # Find values of all indices per assemblage
  fd_ind_values_sharks.CH <- alpha_fd_indices_sharks.CH$"functional_diversity_indices"
  
  ## Form dataframe
  FD.CH <- as.data.frame(fd_ind_values_sharks.CH) %>% 
    tibble::rownames_to_column("Epoch")
  FD.CH$Epoch <- ordered(FD.CH$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                               "Pliocene","Pleistocene","Recent"))
  
  # Form list to merge datasets
  FDind.CH = list(FEmetrics.CH,FD.CH)
  
  # Output
  FDind.CH %>% 
    reduce(inner_join, by = "Epoch")
  
  #close loop
})

stopImplicitCluster()

# Merge lists into 1 dataframe

res_df_CH <- res.CH %>% 
  bind_rows()

# Remove sp_richn, which is the same as nb_sp - redundant
res_df_CH$sp_richn<-NULL

# Format res for FD metrics
Res_FDmetrics_CH<- res_df_CH %>% 
  select(Epoch:fspe)
save(Res_FDmetrics_CH, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/CH_variation_metrics.RData")

# Melt data
FDmetrics_long_CH<- melt(Res_FDmetrics_CH, id.vars= "Epoch")
save(FDmetrics_long_CH, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/CH_variation_long_metrics.RData")

# Form dataframe of mean, median and standard deviation of all FD metrics
Taxon_CH <- Res_FDmetrics_CH %>% 
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

Taxon_CH$Epoch <- ordered(Taxon_CH$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                   "Pliocene","Pleistocene","Recent"))
Taxon_CH$Spp <- Spp_richness

# Save iteration data
save(Taxon_CH, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/CH_Taxon_metrics.RData")


# Analyses without CW
# Form FTU unique combinations, with epoch included
abun_CW <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CE,LC,XO,LO,Epoch,Epoch_earliest,Epoch_latest,Current_status)

# Order taxa alphabetically like before and order Epochs
abun_CW <- abun_CW[order(abun_CW$Taxon_corrected),]
abun_CW$Epoch <- ordered(abun_CW$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                 "Pliocene","Pleistocene", "Recent","nonmatch"))

# Form FTUs column
abun_CW$FTU<-paste(abun_CW$Taxon_corrected, abun_CW$CH, abun_CW$CE, abun_CW$LC, abun_CW$XO, abun_CW$LO, sep="+")

# Cast into wide-format data
wp.abun_CW <- abun_CW %>% 
  pivot_wider(names_from = Epoch, values_from = FTU)

wp.abun_CW <- wp.abun_CW[, c("Taxon_corrected","Epoch_earliest","Epoch_latest","Current_status",
                             "Palaeocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene","Recent","nonmatch")]


# Mark 0 and 1 for absence and presence data in all epochs
wp.abun_CW$Palaeocene[is.na(wp.abun_CW$Palaeocene)] <- "0" 
wp.abun_CW$Eocene[is.na(wp.abun_CW$Eocene)] <- "0"
wp.abun_CW$Oligocene[is.na(wp.abun_CW$Oligocene)] <- "0"
wp.abun_CW$Miocene[is.na(wp.abun_CW$Miocene)] <- "0"
wp.abun_CW$Pliocene[is.na(wp.abun_CW$Pliocene)] <- "0"
wp.abun_CW$Pleistocene[is.na(wp.abun_CW$Pleistocene)] <- "0"
wp.abun_CW$Recent[is.na(wp.abun_CW$Recent)] <- "0"
wp.abun_CW$nonmatch[is.na(wp.abun_CW$nonmatch)] <- "0"

wp.abun_CW <- wp.abun_CW %>% 
  mutate(Palaeocene = replace(Palaeocene,Palaeocene!="0", "1")) %>%
  mutate(Eocene = replace(Eocene,Eocene!="0", "1")) %>% 
  mutate(Oligocene = replace(Oligocene,Oligocene!="0", "1")) %>% 
  mutate(Miocene = replace(Miocene,Miocene!="0", "1")) %>%  
  mutate(Pliocene = replace(Pliocene,Pliocene!="0", "1")) %>%
  mutate(Pleistocene = replace(Pleistocene,Pleistocene!="0", "1")) %>% 
  mutate(Recent = replace(Recent,Recent!="0", "1")) %>%
  mutate(nonmatch = replace(nonmatch,nonmatch!="0", "1"))

## Fill in non-matches
wp.abun_CW <- wp.abun_CW %>% 
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
wp.abun_CW <- wp.abun_CW %>% 
  rowwise() %>% 
  mutate(Eocene = case_when(Palaeocene == "1" && Current_status == "Extant"~"1", TRUE ~ Eocene)) %>% 
  mutate(Oligocene = case_when(Eocene == "1" && Current_status == "Extant"~"1", TRUE ~ Oligocene)) %>% 
  mutate(Miocene = case_when(Oligocene == "1" && Current_status == "Extant"~"1", TRUE ~ Miocene)) %>% 
  mutate(Pliocene = case_when(Miocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pliocene)) %>% 
  mutate(Pleistocene = case_when(Pliocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pleistocene)) %>% 
  mutate(Recent = case_when(Pleistocene == "1" && Current_status == "Extant"~"1", TRUE ~ Recent))

baskets.range <- wp.abun_CW %>% 
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

Spp_richness <- c(length(unique(baskets.epoch$Palaeocene)),
                  length(unique(baskets.epoch$Eocene)),
                  length(unique(baskets.epoch$Oligocene)),
                  length(unique(baskets.epoch$Miocene)),
                  length(unique(baskets.epoch$Pliocene)),
                  length(unique(baskets.epoch$Pleistocene)),
                  length(unique(baskets.epoch$Recent)))

# Form functional taxonomic units (FTUs), grouping by species to account for intraspecific variability
FTU_CW <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CE,LC,XO,LO)

# Order FTUs alphabetically to help ensure matches between species-trait matrix and occurrence matrix later
FTU_CW <- FTU_CW[order(FTU_CW$Taxon_corrected),]

# Filter taxa that are NOT present in "selected_taxa"
FTU_filtered <- FTU_CW %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

## Run iterations & calculate mean by taxon
# Set up model to run on n iterations
registerDoParallel(cores = 10)
getDoParWorkers()
res.CW<-NULL
res.CW<- lapply(1:1000, function(x){
  # Form sliced species-trait matrix
  sharks_traits.CW <- FTU_filtered %>%
    group_by(Taxon_corrected) %>%
    slice_sample(n=1, replace = FALSE) %>% 
    remove_rownames %>% 
    column_to_rownames(var="Taxon_corrected")
  
  # Transpose to form final occurrence matrix
  baskets_sharks_weights.CW <- t(baskets)
  
  # Make matrix and numeric
  baskets_sharks_weights.CW <- data.matrix(baskets_sharks_weights.CW, rownames.force = NA)
  class(baskets_sharks_weights.CW) <- "numeric"
  
  # Summarise & extract occurrences
  asb_sp_sharks_summ.CW <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights.CW)
  
  asb_sp_sharks_occ.CW <- asb_sp_sharks_summ.CW$asb_sp_occ
  
  ############################################################################################
  # Make trait_cat matrix
  sharks_traits_cat.CW <- tibble(trait_name = c("CH", "CE", "LC", "XO", "LO"),
                                 trait_type = c("O", "N", "N", "N", "N"),
                                 trait_weight = c(0.5,1,1,0.33,0.67))
  
  ############################################################################################
  # Form the functional entities
  sp_to_fe_sharks.CW <- mFD::sp.to.fe(
    sp_tr       = sharks_traits.CW, 
    tr_cat      = sharks_traits_cat.CW, 
    fe_nm_type  = "fe_rank", 
    check_input = TRUE)
  
  # Calculate and display metrics
  alpha_fd_fe_sharks.CW <- alpha.fd.fe(asb_sp_sharks_occ.CW,
                                       sp_to_fe_sharks.CW,
                                       ind_nm = c("fred", "fored", "fvuln"),
                                       check_input = TRUE,
                                       details_returned = TRUE)
  
  # Form dataframe
  FEmetrics.CW <- as.data.frame(alpha_fd_fe_sharks.CW$asb_fdfe) %>% 
    tibble::rownames_to_column("Epoch")
  FEmetrics.CW$Epoch <- ordered(FEmetrics.CW$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                             "Pliocene","Pleistocene","Recent"))
  
  # Construct trait distance matrix
  sp_dist_sharks.CW <- mFD::funct.dist(
    sp_tr         = sharks_traits.CW,
    tr_cat        = sharks_traits_cat.CW,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  # Assess quality
  fspaces_quality_sharks.CW <- mFD::quality.fspaces(
    sp_dist             = sp_dist_sharks.CW,
    maxdim_pcoa         = 10,
    deviation_weighting = c("absolute","squared"),
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  # Return coordinates of each PC/axis
  sp_faxes_coord_sharks.CW <- fspaces_quality_sharks.CW$"details_fspaces"$"sp_pc_coord"
  
  # Form alpha diversity indices
  alpha_fd_indices_sharks.CW <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_faxes_coord_sharks.CW[ , c("PC1", "PC2", "PC3")],
    asb_sp_w         = baskets_sharks_weights.CW,
    ind_vect         = c("fric", "fori", "fspe"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  # Find values of all indices per assemblage
  fd_ind_values_sharks.CW <- alpha_fd_indices_sharks.CW$"functional_diversity_indices"
  
  ## Form dataframe
  FD.CW <- as.data.frame(fd_ind_values_sharks.CW) %>% 
    tibble::rownames_to_column("Epoch")
  FD.CW$Epoch <- ordered(FD.CW$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                               "Pliocene","Pleistocene","Recent"))
  
  # Form list to merge datasets
  FDind.CW = list(FEmetrics.CW,FD.CW)
  
  # Output
  FDind.CW %>% 
    reduce(inner_join, by = "Epoch")
  
  #close loop
})

stopImplicitCluster()

# Merge lists into 1 dataframe
res_df_CW <- res.CW %>% 
  bind_rows()

# Remove sp_richn, which is the same as nb_sp - redundant
res_df_CW$sp_richn<-NULL

# Format res for FD metrics
Res_FDmetrics_CW<- res_df_CW %>% 
  select(Epoch:fspe)
save(Res_FDmetrics_CW, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/CW_variation_metrics.RData")

# Melt data
FDmetrics_long_CW<- melt(Res_FDmetrics_CW, id.vars= "Epoch")
save(FDmetrics_long_CW, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/CW_variation_long_metrics.RData")

# Form dataframe of mean, median and standard deviation of all FD metrics
Taxon_CW <- Res_FDmetrics_CW %>% 
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

Taxon_CW$Epoch <- ordered(Taxon_CW$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                   "Pliocene","Pleistocene","Recent"))
Taxon_CW$Spp <- Spp_richness

# Save iteration data
save(Taxon_CW, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/CW_Taxon_metrics.RData")

# Analyses without CE
# Form FTU unique combinations, with epoch included
abun_CE <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,LC,XO,LO,Epoch,Epoch_earliest,Epoch_latest,Current_status)

# Order taxa alphabetically like before and order Epochs
abun_CE <- abun_CE[order(abun_CE$Taxon_corrected),]
abun_CE$Epoch <- ordered(abun_CE$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                 "Pliocene","Pleistocene", "Recent","nonmatch"))

# Form FTUs column
abun_CE$FTU<-paste(abun_CE$Taxon_corrected, abun_CE$CH, abun_CE$CW, abun_CE$LC, abun_CE$XO, abun_CE$LO, sep="+")

# Cast into wide-format data
wp.abun_CE <- abun_CE %>% 
  pivot_wider(names_from = Epoch, values_from = FTU)

wp.abun_CE <- wp.abun_CE[, c("Taxon_corrected","Epoch_earliest","Epoch_latest","Current_status",
                             "Palaeocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene","Recent","nonmatch")]


# Mark 0 and 1 for absence and presence data in all epochs
wp.abun_CE$Palaeocene[is.na(wp.abun_CE$Palaeocene)] <- "0" 
wp.abun_CE$Eocene[is.na(wp.abun_CE$Eocene)] <- "0"
wp.abun_CE$Oligocene[is.na(wp.abun_CE$Oligocene)] <- "0"
wp.abun_CE$Miocene[is.na(wp.abun_CE$Miocene)] <- "0"
wp.abun_CE$Pliocene[is.na(wp.abun_CE$Pliocene)] <- "0"
wp.abun_CE$Pleistocene[is.na(wp.abun_CE$Pleistocene)] <- "0"
wp.abun_CE$Recent[is.na(wp.abun_CE$Recent)] <- "0"
wp.abun_CE$nonmatch[is.na(wp.abun_CE$nonmatch)] <- "0"

wp.abun_CE <- wp.abun_CE %>% 
  mutate(Palaeocene = replace(Palaeocene,Palaeocene!="0", "1")) %>%
  mutate(Eocene = replace(Eocene,Eocene!="0", "1")) %>% 
  mutate(Oligocene = replace(Oligocene,Oligocene!="0", "1")) %>% 
  mutate(Miocene = replace(Miocene,Miocene!="0", "1")) %>%  
  mutate(Pliocene = replace(Pliocene,Pliocene!="0", "1")) %>%
  mutate(Pleistocene = replace(Pleistocene,Pleistocene!="0", "1")) %>% 
  mutate(Recent = replace(Recent,Recent!="0", "1")) %>%
  mutate(nonmatch = replace(nonmatch,nonmatch!="0", "1"))

## Fill in non-matches
wp.abun_CE <- wp.abun_CE %>% 
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
wp.abun_CE <- wp.abun_CE %>% 
  rowwise() %>% 
  mutate(Eocene = case_when(Palaeocene == "1" && Current_status == "Extant"~"1", TRUE ~ Eocene)) %>% 
  mutate(Oligocene = case_when(Eocene == "1" && Current_status == "Extant"~"1", TRUE ~ Oligocene)) %>% 
  mutate(Miocene = case_when(Oligocene == "1" && Current_status == "Extant"~"1", TRUE ~ Miocene)) %>% 
  mutate(Pliocene = case_when(Miocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pliocene)) %>% 
  mutate(Pleistocene = case_when(Pliocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pleistocene)) %>% 
  mutate(Recent = case_when(Pleistocene == "1" && Current_status == "Extant"~"1", TRUE ~ Recent))

baskets.range <- wp.abun_CE %>% 
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

Spp_richness <- c(length(unique(baskets.epoch$Palaeocene)),
                  length(unique(baskets.epoch$Eocene)),
                  length(unique(baskets.epoch$Oligocene)),
                  length(unique(baskets.epoch$Miocene)),
                  length(unique(baskets.epoch$Pliocene)),
                  length(unique(baskets.epoch$Pleistocene)),
                  length(unique(baskets.epoch$Recent)))

# Form functional taxonomic units (FTUs), grouping by species to account for intraspecific variability
FTU_CE <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,LC,XO,LO)

# Order FTUs alphabetically to help ensure matches between species-trait matrix and occurrence matrix later
FTU_CE <- FTU_CE[order(FTU_CE$Taxon_corrected),]

# Filter taxa that are NOT present in "selected_taxa"
FTU_filtered <- FTU_CE %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

## Run iterations & calculate mean by taxon
# Set up model to run on n iterations
registerDoParallel(cores = 10)
getDoParWorkers()
res.CE<-NULL
res.CE<- lapply(1:1000, function(x){
  # Form sliced species-trait matrix
  sharks_traits.CE <- FTU_filtered %>%
    group_by(Taxon_corrected) %>%
    slice_sample(n=1, replace = FALSE) %>% 
    remove_rownames %>% 
    column_to_rownames(var="Taxon_corrected")
  
  # Transpose to form final occurrence matrix
  baskets_sharks_weights.CE <- t(baskets)
  
  # Make matrix and numeric
  baskets_sharks_weights.CE <- data.matrix(baskets_sharks_weights.CE, rownames.force = NA)
  class(baskets_sharks_weights.CE) <- "numeric"
  
  # Summarise & extract occurrences
  asb_sp_sharks_summ.CE <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights.CE)
  
  asb_sp_sharks_occ.CE <- asb_sp_sharks_summ.CE$asb_sp_occ
  
  ############################################################################################
  # Make trait_cat matrix
  sharks_traits_cat.CE <- tibble(trait_name = c("CH", "CW", "LC", "XO", "LO"),
                                 trait_type = c("O", "O", "N", "N", "N"),
                                 trait_weight = c(0.5,0.5,1,0.33,0.67))
  
  ############################################################################################
  # Form the functional entities
  sp_to_fe_sharks.CE <- mFD::sp.to.fe(
    sp_tr       = sharks_traits.CE, 
    tr_cat      = sharks_traits_cat.CE, 
    fe_nm_type  = "fe_rank", 
    check_input = TRUE)
  
  # Calculate and display metrics
  alpha_fd_fe_sharks.CE <- alpha.fd.fe(asb_sp_sharks_occ.CE,
                                       sp_to_fe_sharks.CE,
                                       ind_nm = c("fred", "fored", "fvuln"),
                                       check_input = TRUE,
                                       details_returned = TRUE)
  
  # Form dataframe
  FEmetrics.CE <- as.data.frame(alpha_fd_fe_sharks.CE$asb_fdfe) %>% 
    tibble::rownames_to_column("Epoch")
  FEmetrics.CE$Epoch <- ordered(FEmetrics.CE$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                             "Pliocene","Pleistocene","Recent"))
  
  # Construct trait distance matrix
  sp_dist_sharks.CE <- mFD::funct.dist(
    sp_tr         = sharks_traits.CE,
    tr_cat        = sharks_traits_cat.CE,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  # Assess quality
  fspaces_quality_sharks.CE <- mFD::quality.fspaces(
    sp_dist             = sp_dist_sharks.CE,
    maxdim_pcoa         = 10,
    deviation_weighting = c("absolute","squared"),
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  # Return coordinates of each PC/axis
  sp_faxes_coord_sharks.CE <- fspaces_quality_sharks.CE$"details_fspaces"$"sp_pc_coord"
  
  # Form alpha diversity indices
  alpha_fd_indices_sharks.CE <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_faxes_coord_sharks.CE[ , c("PC1", "PC2", "PC3")],
    asb_sp_w         = baskets_sharks_weights.CE,
    ind_vect         = c("fric", "fori", "fspe"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  # Find values of all indices per assemblage
  fd_ind_values_sharks.CE <- alpha_fd_indices_sharks.CE$"functional_diversity_indices"
  
  ## Form dataframe
  FD.CE <- as.data.frame(fd_ind_values_sharks.CE) %>% 
    tibble::rownames_to_column("Epoch")
  FD.CE$Epoch <- ordered(FD.CE$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                               "Pliocene","Pleistocene","Recent"))
  
  # Form list to merge datasets
  FDind.CE = list(FEmetrics.CE,FD.CE)
  
  # Output
  FDind.CE %>% 
    reduce(inner_join, by = "Epoch")
  
  #close loop
})

stopImplicitCluster()
#Merge lists into 1 dataframe
res_df_CE <- res.CE %>% 
  bind_rows()

# Remove sp_richn, which is the same as nb_sp - redundant
res_df_CE$sp_richn<-NULL

# Format res for FD metrics
Res_FDmetrics_CE<- res_df_CE %>% 
  select(Epoch:fspe)
save(Res_FDmetrics_CE, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/CE_variation_metrics.RData")

# Melt data
FDmetrics_long_CE<- melt(Res_FDmetrics_CE, id.vars= "Epoch")
save(FDmetrics_long_CE, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/CE_variation_long_metrics.RData")

# Form dataframe of mean, median and standard deviation of all FD metrics
Taxon_CE <- Res_FDmetrics_CE %>% 
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

Taxon_CE$Epoch <- ordered(Taxon_CW$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                   "Pliocene","Pleistocene","Recent"))
Taxon_CE$Spp <- Spp_richness

# Save iteration data
save(Taxon_CE, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/CE_Taxon_metrics.RData")

# Analyses without LC
# Form FTU unique combinations, with epoch included
abun_LC <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,CE,XO,LO,Epoch,Epoch_earliest,Epoch_latest,Current_status)

# Order taxa alphabetically like before and order Epochs
abun_LC <- abun_LC[order(abun_LC$Taxon_corrected),]
abun_LC$Epoch <- ordered(abun_LC$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                 "Pliocene","Pleistocene", "Recent","nonmatch"))

# Form FTUs column
abun_LC$FTU<-paste(abun_LC$Taxon_corrected, abun_LC$CH, abun_LC$CW, abun_LC$CE, abun_LC$XO, abun_LC$LO, sep="+")

# Cast into wide-format data
wp.abun_LC <- abun_LC %>% 
  pivot_wider(names_from = Epoch, values_from = FTU)

wp.abun_LC <- wp.abun_LC[, c("Taxon_corrected","Epoch_earliest","Epoch_latest","Current_status",
                             "Palaeocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene","Recent","nonmatch")]


# Mark 0 and 1 for absence and presence data in all epochs
wp.abun_LC$Palaeocene[is.na(wp.abun_LC$Palaeocene)] <- "0" 
wp.abun_LC$Eocene[is.na(wp.abun_LC$Eocene)] <- "0"
wp.abun_LC$Oligocene[is.na(wp.abun_LC$Oligocene)] <- "0"
wp.abun_LC$Miocene[is.na(wp.abun_LC$Miocene)] <- "0"
wp.abun_LC$Pliocene[is.na(wp.abun_LC$Pliocene)] <- "0"
wp.abun_LC$Pleistocene[is.na(wp.abun_LC$Pleistocene)] <- "0"
wp.abun_LC$Recent[is.na(wp.abun_LC$Recent)] <- "0"
wp.abun_LC$nonmatch[is.na(wp.abun_LC$nonmatch)] <- "0"

wp.abun_LC <- wp.abun_LC %>% 
  mutate(Palaeocene = replace(Palaeocene,Palaeocene!="0", "1")) %>%
  mutate(Eocene = replace(Eocene,Eocene!="0", "1")) %>% 
  mutate(Oligocene = replace(Oligocene,Oligocene!="0", "1")) %>% 
  mutate(Miocene = replace(Miocene,Miocene!="0", "1")) %>%  
  mutate(Pliocene = replace(Pliocene,Pliocene!="0", "1")) %>%
  mutate(Pleistocene = replace(Pleistocene,Pleistocene!="0", "1")) %>% 
  mutate(Recent = replace(Recent,Recent!="0", "1")) %>%
  mutate(nonmatch = replace(nonmatch,nonmatch!="0", "1"))

## Fill in non-matches
wp.abun_LC <- wp.abun_LC %>% 
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
wp.abun_LC <- wp.abun_LC %>% 
  rowwise() %>% 
  mutate(Eocene = case_when(Palaeocene == "1" && Current_status == "Extant"~"1", TRUE ~ Eocene)) %>% 
  mutate(Oligocene = case_when(Eocene == "1" && Current_status == "Extant"~"1", TRUE ~ Oligocene)) %>% 
  mutate(Miocene = case_when(Oligocene == "1" && Current_status == "Extant"~"1", TRUE ~ Miocene)) %>% 
  mutate(Pliocene = case_when(Miocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pliocene)) %>% 
  mutate(Pleistocene = case_when(Pliocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pleistocene)) %>% 
  mutate(Recent = case_when(Pleistocene == "1" && Current_status == "Extant"~"1", TRUE ~ Recent))

baskets.range <- wp.abun_LC %>% 
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

Spp_richness <- c(length(unique(baskets.epoch$Palaeocene)),
                  length(unique(baskets.epoch$Eocene)),
                  length(unique(baskets.epoch$Oligocene)),
                  length(unique(baskets.epoch$Miocene)),
                  length(unique(baskets.epoch$Pliocene)),
                  length(unique(baskets.epoch$Pleistocene)),
                  length(unique(baskets.epoch$Recent)))

# Form functional taxonomic units (FTUs), grouping by species to account for intraspecific variability
FTU_LC <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,CE,XO,LO)

# Order FTUs alphabetically to help ensure matches between species-trait matrix and occurrence matrix later
FTU_LC <- FTU_LC[order(FTU_LC$Taxon_corrected),]

# Filter taxa that are NOT present in "selected_taxa"
FTU_filtered <- FTU_LC %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

## Run iterations & calculate mean by taxon
# Set up model to run on n iterations
registerDoParallel(cores = 10)
getDoParWorkers()
res.LC<-NULL
res.LC<- lapply(1:1000, function(x){
  # Form sliced species-trait matrix
  sharks_traits.LC <- FTU_filtered %>%
    group_by(Taxon_corrected) %>%
    slice_sample(n=1, replace = FALSE) %>% 
    remove_rownames %>% 
    column_to_rownames(var="Taxon_corrected")
  
  # Transpose to form final occurrence matrix
  baskets_sharks_weights.LC <- t(baskets)
  
  # Make matrix and numeric
  baskets_sharks_weights.LC <- data.matrix(baskets_sharks_weights.LC, rownames.force = NA)
  class(baskets_sharks_weights.LC) <- "numeric"
  
  # Summarise & extract occurrences
  asb_sp_sharks_summ.LC <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights.LC)
  
  asb_sp_sharks_occ.LC <- asb_sp_sharks_summ.LC$asb_sp_occ
  
  ############################################################################################
  # Make trait_cat matrix
  sharks_traits_cat.LC <- tibble(trait_name = c("CH", "CW", "CE", "XO", "LO"),
                                 trait_type = c("O", "O", "N", "N", "N"),
                                 trait_weight = c(0.5,0.5,1,0.33,0.67))
  
  ############################################################################################
  # Form the functional entities
  sp_to_fe_sharks.LC <- mFD::sp.to.fe(
    sp_tr       = sharks_traits.LC, 
    tr_cat      = sharks_traits_cat.LC, 
    fe_nm_type  = "fe_rank", 
    check_input = TRUE)
  
  # Calculate and display metrics
  alpha_fd_fe_sharks.LC <- alpha.fd.fe(asb_sp_sharks_occ.LC,
                                       sp_to_fe_sharks.LC,
                                       ind_nm = c("fred", "fored", "fvuln"),
                                       check_input = TRUE,
                                       details_returned = TRUE)
  
  # Form dataframe
  FEmetrics.LC <- as.data.frame(alpha_fd_fe_sharks.LC$asb_fdfe) %>% 
    tibble::rownames_to_column("Epoch")
  FEmetrics.LC$Epoch <- ordered(FEmetrics.LC$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                             "Pliocene","Pleistocene","Recent"))
  
  # Construct trait distance matrix
  sp_dist_sharks.LC <- mFD::funct.dist(
    sp_tr         = sharks_traits.LC,
    tr_cat        = sharks_traits_cat.LC,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  # Assess quality
  fspaces_quality_sharks.LC <- mFD::quality.fspaces(
    sp_dist             = sp_dist_sharks.LC,
    maxdim_pcoa         = 10,
    deviation_weighting = c("absolute","squared"),
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  # Return coordinates of each PC/axis
  sp_faxes_coord_sharks.LC <- fspaces_quality_sharks.LC$"details_fspaces"$"sp_pc_coord"
  
  # Form alpha diversity indices
  alpha_fd_indices_sharks.LC <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_faxes_coord_sharks.LC[ , c("PC1", "PC2", "PC3")],
    asb_sp_w         = baskets_sharks_weights.LC,
    ind_vect         = c("fric", "fori", "fspe"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  # Find values of all indices per assemblage
  fd_ind_values_sharks.LC <- alpha_fd_indices_sharks.LC$"functional_diversity_indices"
  
  ## Form dataframe
  FD.LC <- as.data.frame(fd_ind_values_sharks.LC) %>% 
    tibble::rownames_to_column("Epoch")
  FD.LC$Epoch <- ordered(FD.LC$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                               "Pliocene","Pleistocene","Recent"))
  
  # Form list to merge datasets
  FDind.LC = list(FEmetrics.LC,FD.LC)
  
  # Output
  FDind.LC %>% 
    reduce(inner_join, by = "Epoch")
  
  #close loop
})

stopImplicitCluster()
#Merge lists into 1 dataframe
res_df_LC <- res.LC %>% 
  bind_rows()

# Remove sp_richn, which is the same as nb_sp - redundant
res_df_LC$sp_richn<-NULL

# Format res for FD metrics
Res_FDmetrics_LC<- res_df_LC %>% 
  select(Epoch:fspe)
save(Res_FDmetrics_LC, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/LC_variation_metrics.RData")

# Melt data
FDmetrics_long_LC<- melt(Res_FDmetrics_LC, id.vars= "Epoch")
save(FDmetrics_long_LC, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/LC_variation_long_metrics.RData")

# Form dataframe of mean, median and standard deviation of all FD metrics
Taxon_LC <- Res_FDmetrics_LC %>% 
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

Taxon_LC$Epoch <- ordered(Taxon_CW$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                   "Pliocene","Pleistocene","Recent"))
Taxon_LC$Spp <- Spp_richness

# Save iteration data
save(Taxon_LC, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/LC_Taxon_metrics.RData")

# Analyses without XO
# Form FTU unique combinations, with epoch included
abun_XO <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,CE,LC,LO,Epoch,Epoch_earliest,Epoch_latest,Current_status)

# Order taxa alphabetically like before and order Epochs
abun_XO <- abun_XO[order(abun_XO$Taxon_corrected),]
abun_XO$Epoch <- ordered(abun_XO$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                 "Pliocene","Pleistocene", "Recent","nonmatch"))

# Form FTUs column
abun_XO$FTU<-paste(abun_XO$Taxon_corrected, abun_XO$CH, abun_XO$CW, abun_XO$CE, abun_XO$LC, abun_XO$LO, sep="+")

# Cast into wide-format data
wp.abun_XO <- abun_XO %>% 
  pivot_wider(names_from = Epoch, values_from = FTU)

wp.abun_XO <- wp.abun_XO[, c("Taxon_corrected","Epoch_earliest","Epoch_latest","Current_status",
                             "Palaeocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene","Recent","nonmatch")]


# Mark 0 and 1 for absence and presence data in all epochs
wp.abun_XO$Palaeocene[is.na(wp.abun_XO$Palaeocene)] <- "0" 
wp.abun_XO$Eocene[is.na(wp.abun_XO$Eocene)] <- "0"
wp.abun_XO$Oligocene[is.na(wp.abun_XO$Oligocene)] <- "0"
wp.abun_XO$Miocene[is.na(wp.abun_XO$Miocene)] <- "0"
wp.abun_XO$Pliocene[is.na(wp.abun_XO$Pliocene)] <- "0"
wp.abun_XO$Pleistocene[is.na(wp.abun_XO$Pleistocene)] <- "0"
wp.abun_XO$Recent[is.na(wp.abun_XO$Recent)] <- "0"
wp.abun_XO$nonmatch[is.na(wp.abun_XO$nonmatch)] <- "0"

wp.abun_XO <- wp.abun_XO %>% 
  mutate(Palaeocene = replace(Palaeocene,Palaeocene!="0", "1")) %>%
  mutate(Eocene = replace(Eocene,Eocene!="0", "1")) %>% 
  mutate(Oligocene = replace(Oligocene,Oligocene!="0", "1")) %>% 
  mutate(Miocene = replace(Miocene,Miocene!="0", "1")) %>%  
  mutate(Pliocene = replace(Pliocene,Pliocene!="0", "1")) %>%
  mutate(Pleistocene = replace(Pleistocene,Pleistocene!="0", "1")) %>% 
  mutate(Recent = replace(Recent,Recent!="0", "1")) %>%
  mutate(nonmatch = replace(nonmatch,nonmatch!="0", "1"))

## Fill in non-matches
wp.abun_XO <- wp.abun_XO %>% 
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
wp.abun_XO <- wp.abun_XO %>% 
  rowwise() %>% 
  mutate(Eocene = case_when(Palaeocene == "1" && Current_status == "Extant"~"1", TRUE ~ Eocene)) %>% 
  mutate(Oligocene = case_when(Eocene == "1" && Current_status == "Extant"~"1", TRUE ~ Oligocene)) %>% 
  mutate(Miocene = case_when(Oligocene == "1" && Current_status == "Extant"~"1", TRUE ~ Miocene)) %>% 
  mutate(Pliocene = case_when(Miocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pliocene)) %>% 
  mutate(Pleistocene = case_when(Pliocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pleistocene)) %>% 
  mutate(Recent = case_when(Pleistocene == "1" && Current_status == "Extant"~"1", TRUE ~ Recent))

baskets.range <- wp.abun_XO %>% 
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

Spp_richness <- c(length(unique(baskets.epoch$Palaeocene)),
                  length(unique(baskets.epoch$Eocene)),
                  length(unique(baskets.epoch$Oligocene)),
                  length(unique(baskets.epoch$Miocene)),
                  length(unique(baskets.epoch$Pliocene)),
                  length(unique(baskets.epoch$Pleistocene)),
                  length(unique(baskets.epoch$Recent)))

# Form functional taxonomic units (FTUs), grouping by species to account for intraspecific variability
FTU_XO <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,CE,LC,LO)

# Order FTUs alphabetically to help ensure matches between species-trait matrix and occurrence matrix later
FTU_XO <- FTU_XO[order(FTU_XO$Taxon_corrected),]

# Filter taxa that are NOT present in "selected_taxa"
FTU_filtered <- FTU_XO %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

## Run iterations & calculate mean by taxon
# Set up model to run on n iterations
registerDoParallel(cores = 10)
getDoParWorkers()
res.XO<-NULL
res.XO<- lapply(1:1000, function(x){
  # Form sliced species-trait matrix
  sharks_traits.XO <- FTU_filtered %>%
    group_by(Taxon_corrected) %>%
    slice_sample(n=1, replace = FALSE) %>% 
    remove_rownames %>% 
    column_to_rownames(var="Taxon_corrected")
  
  # Transpose to form final occurrence matrix
  baskets_sharks_weights.XO <- t(baskets)
  
  # Make matrix and numeric
  baskets_sharks_weights.XO <- data.matrix(baskets_sharks_weights.XO, rownames.force = NA)
  class(baskets_sharks_weights.XO) <- "numeric"
  
  # Summarise & extract occurrences
  asb_sp_sharks_summ.XO <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights.XO)
  
  asb_sp_sharks_occ.XO <- asb_sp_sharks_summ.XO$asb_sp_occ
  
  ############################################################################################
  # Make trait_cat matrix
  sharks_traits_cat.XO <- tibble(trait_name = c("CH", "CW", "CE", "LC", "LO"),
                                 trait_type = c("O", "O", "N", "N", "N"),
                                 trait_weight = c(0.5,0.5,1,1,0.67))
  
  ############################################################################################
  # Form the functional entities
  sp_to_fe_sharks.XO <- mFD::sp.to.fe(
    sp_tr       = sharks_traits.XO, 
    tr_cat      = sharks_traits_cat.XO, 
    fe_nm_type  = "fe_rank", 
    check_input = TRUE)
  
  # Calculate and display metrics
  alpha_fd_fe_sharks.XO <- alpha.fd.fe(asb_sp_sharks_occ.XO,
                                       sp_to_fe_sharks.XO,
                                       ind_nm = c("fred", "fored", "fvuln"),
                                       check_input = TRUE,
                                       details_returned = TRUE)
  
  # Form dataframe
  FEmetrics.XO <- as.data.frame(alpha_fd_fe_sharks.XO$asb_fdfe) %>% 
    tibble::rownames_to_column("Epoch")
  FEmetrics.XO$Epoch <- ordered(FEmetrics.XO$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                             "Pliocene","Pleistocene","Recent"))
  
  # Construct trait distance matrix
  sp_dist_sharks.XO <- mFD::funct.dist(
    sp_tr         = sharks_traits.XO,
    tr_cat        = sharks_traits_cat.XO,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  # Assess quality
  fspaces_quality_sharks.XO <- mFD::quality.fspaces(
    sp_dist             = sp_dist_sharks.XO,
    maxdim_pcoa         = 10,
    deviation_weighting = c("absolute","squared"),
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  # Return coordinates of each PC/axis
  sp_faxes_coord_sharks.XO <- fspaces_quality_sharks.XO$"details_fspaces"$"sp_pc_coord"
  
  # Form alpha diversity indices
  alpha_fd_indices_sharks.XO <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_faxes_coord_sharks.XO[ , c("PC1", "PC2", "PC3")],
    asb_sp_w         = baskets_sharks_weights.XO,
    ind_vect         = c("fric", "fori", "fspe"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  # Find values of all indices per assemblage
  fd_ind_values_sharks.XO <- alpha_fd_indices_sharks.XO$"functional_diversity_indices"
  
  ## Form dataframe
  FD.XO <- as.data.frame(fd_ind_values_sharks.XO) %>% 
    tibble::rownames_to_column("Epoch")
  FD.XO$Epoch <- ordered(FD.XO$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                               "Pliocene","Pleistocene","Recent"))
  
  # Form list to merge datasets
  FDind.XO = list(FEmetrics.XO,FD.XO)
  
  # Output
  FDind.XO %>% 
    reduce(inner_join, by = "Epoch")
  
  #close loop
})

stopImplicitCluster()

# Merge into 1 dataframe
res_df_XO <- res.XO %>% 
  bind_rows()

# Remove sp_richn, which is the same as nb_sp - redundant
res_df_XO$sp_richn<-NULL

# Format res for FD metrics
Res_FDmetrics_XO<- res_df_XO %>% 
  select(Epoch:fspe)
save(Res_FDmetrics_XO, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/XO_variation_metrics.RData")

# Melt data
FDmetrics_long_XO<- melt(Res_FDmetrics_XO, id.vars= "Epoch")
save(FDmetrics_long_XO, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/XO_variation_long_metrics.RData")

# Form dataframe of mean, median and standard deviation of all FD metrics
Taxon_XO <- Res_FDmetrics_XO %>% 
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

Taxon_XO$Epoch <- ordered(Taxon_CW$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                   "Pliocene","Pleistocene","Recent"))
Taxon_XO$Spp <- Spp_richness

# Save iteration data
save(Taxon_XO, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/XO_Taxon_metrics.RData")

# Analyses without LO
# Form FTU unique combinations, with epoch included
abun_LO <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,CE,LC,XO,Epoch,Epoch_earliest,Epoch_latest,Current_status)

# Order taxa alphabetically like before and order Epochs
abun_LO <- abun_LO[order(abun_LO$Taxon_corrected),]
abun_LO$Epoch <- ordered(abun_LO$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                 "Pliocene","Pleistocene", "Recent","nonmatch"))

# Form FTUs column
abun_LO$FTU<-paste(abun_LO$Taxon_corrected, abun_LO$CH, abun_LO$CW, abun_LO$CE, abun_LO$LC, abun_LO$XO, sep="+")

# Cast into wide-format data
wp.abun_LO <- abun_LO %>% 
  pivot_wider(names_from = Epoch, values_from = FTU)

wp.abun_LO <- wp.abun_LO[, c("Taxon_corrected","Epoch_earliest","Epoch_latest","Current_status",
                             "Palaeocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene","Recent","nonmatch")]


# Mark 0 and 1 for absence and presence data in all epochs
wp.abun_LO$Palaeocene[is.na(wp.abun_LO$Palaeocene)] <- "0" 
wp.abun_LO$Eocene[is.na(wp.abun_LO$Eocene)] <- "0"
wp.abun_LO$Oligocene[is.na(wp.abun_LO$Oligocene)] <- "0"
wp.abun_LO$Miocene[is.na(wp.abun_LO$Miocene)] <- "0"
wp.abun_LO$Pliocene[is.na(wp.abun_LO$Pliocene)] <- "0"
wp.abun_LO$Pleistocene[is.na(wp.abun_LO$Pleistocene)] <- "0"
wp.abun_LO$Recent[is.na(wp.abun_LO$Recent)] <- "0"
wp.abun_LO$nonmatch[is.na(wp.abun_LO$nonmatch)] <- "0"

wp.abun_LO <- wp.abun_LO %>% 
  mutate(Palaeocene = replace(Palaeocene,Palaeocene!="0", "1")) %>%
  mutate(Eocene = replace(Eocene,Eocene!="0", "1")) %>% 
  mutate(Oligocene = replace(Oligocene,Oligocene!="0", "1")) %>% 
  mutate(Miocene = replace(Miocene,Miocene!="0", "1")) %>%  
  mutate(Pliocene = replace(Pliocene,Pliocene!="0", "1")) %>%
  mutate(Pleistocene = replace(Pleistocene,Pleistocene!="0", "1")) %>% 
  mutate(Recent = replace(Recent,Recent!="0", "1")) %>%
  mutate(nonmatch = replace(nonmatch,nonmatch!="0", "1"))

## Fill in non-matches
wp.abun_LO <- wp.abun_LO %>% 
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
wp.abun_LO <- wp.abun_LO %>% 
  rowwise() %>% 
  mutate(Eocene = case_when(Palaeocene == "1" && Current_status == "Extant"~"1", TRUE ~ Eocene)) %>% 
  mutate(Oligocene = case_when(Eocene == "1" && Current_status == "Extant"~"1", TRUE ~ Oligocene)) %>% 
  mutate(Miocene = case_when(Oligocene == "1" && Current_status == "Extant"~"1", TRUE ~ Miocene)) %>% 
  mutate(Pliocene = case_when(Miocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pliocene)) %>% 
  mutate(Pleistocene = case_when(Pliocene == "1" && Current_status == "Extant"~"1", TRUE ~ Pleistocene)) %>% 
  mutate(Recent = case_when(Pleistocene == "1" && Current_status == "Extant"~"1", TRUE ~ Recent))

baskets.range <- wp.abun_LO %>% 
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

Spp_richness <- c(length(unique(baskets.epoch$Palaeocene)),
                  length(unique(baskets.epoch$Eocene)),
                  length(unique(baskets.epoch$Oligocene)),
                  length(unique(baskets.epoch$Miocene)),
                  length(unique(baskets.epoch$Pliocene)),
                  length(unique(baskets.epoch$Pleistocene)),
                  length(unique(baskets.epoch$Recent)))

# Form functional taxonomic units (FTUs), grouping by species to account for intraspecific variability
FTU_LO <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,CE,LC,XO)

# Order FTUs alphabetically to help ensure matches between species-trait matrix and occurrence matrix later
FTU_LO <- FTU_LO[order(FTU_LO$Taxon_corrected),]

# Filter taxa that are NOT present in "selected_taxa"
FTU_filtered <- FTU_LO %>%
  anti_join(selected_taxa, by = "Taxon_corrected")

## Run iterations & calculate mean by taxon
# Set up model to run on n iterations
registerDoParallel(cores = 10)
getDoParWorkers()
res.LO<-NULL
res.LO<- lapply(1:1000, function(x){
  # Form sliced species-trait matrix
  sharks_traits.LO <- FTU_filtered %>%
    group_by(Taxon_corrected) %>%
    slice_sample(n=1, replace = FALSE) %>% 
    remove_rownames %>% 
    column_to_rownames(var="Taxon_corrected")
  
  # Transpose to form final occurrence matrix
  baskets_sharks_weights.LO <- t(baskets)
  
  # Make matrix and numeric
  baskets_sharks_weights.LO <- data.matrix(baskets_sharks_weights.LO, rownames.force = NA)
  class(baskets_sharks_weights.LO) <- "numeric"
  
  # Summarise & extract occurrences
  asb_sp_sharks_summ.LO <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights.LO)
  
  asb_sp_sharks_occ.LO <- asb_sp_sharks_summ.LO$asb_sp_occ
  
  ############################################################################################
  # Make trait_cat matrix
  sharks_traits_cat.LO <- tibble(trait_name = c("CH", "CW", "CE", "LC", "XO"),
                                 trait_type = c("O", "O", "N", "N", "N"),
                                 trait_weight = c(0.5,0.5,1,1,0.33))
  
  ############################################################################################
  # Form the functional entities
  sp_to_fe_sharks.LO <- mFD::sp.to.fe(
    sp_tr       = sharks_traits.LO, 
    tr_cat      = sharks_traits_cat.LO, 
    fe_nm_type  = "fe_rank", 
    check_input = TRUE)
  
  # Calculate and display metrics
  alpha_fd_fe_sharks.LO <- alpha.fd.fe(asb_sp_sharks_occ.LO,
                                       sp_to_fe_sharks.LO,
                                       ind_nm = c("fred", "fored", "fvuln"),
                                       check_input = TRUE,
                                       details_returned = TRUE)
  
  # Form dataframe
  FEmetrics.LO <- as.data.frame(alpha_fd_fe_sharks.LO$asb_fdfe) %>% 
    tibble::rownames_to_column("Epoch")
  FEmetrics.LO$Epoch <- ordered(FEmetrics.LO$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                             "Pliocene","Pleistocene","Recent"))
  
  # Construct trait distance matrix
  sp_dist_sharks.LO <- mFD::funct.dist(
    sp_tr         = sharks_traits.LO,
    tr_cat        = sharks_traits_cat.LO,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  # Assess quality
  fspaces_quality_sharks.LO <- mFD::quality.fspaces(
    sp_dist             = sp_dist_sharks.LO,
    maxdim_pcoa         = 10,
    deviation_weighting = c("absolute","squared"),
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  # Return coordinates of each PC/axis
  sp_faxes_coord_sharks.LO <- fspaces_quality_sharks.LO$"details_fspaces"$"sp_pc_coord"
  
  # Form alpha diversity indices
  alpha_fd_indices_sharks.LO <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_faxes_coord_sharks.LO[ , c("PC1", "PC2", "PC3")],
    asb_sp_w         = baskets_sharks_weights.LO,
    ind_vect         = c("fric", "fori", "fspe"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  # Find values of all indices per assemblage
  fd_ind_values_sharks.LO <- alpha_fd_indices_sharks.LO$"functional_diversity_indices"
  
  ## Form dataframe
  FD.LO <- as.data.frame(fd_ind_values_sharks.LO) %>% 
    tibble::rownames_to_column("Epoch")
  FD.LO$Epoch <- ordered(FD.LO$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                               "Pliocene","Pleistocene","Recent"))
  
  # Form list to merge datasets
  FDind.LO = list(FEmetrics.LO,FD.LO)
  
  # Output
  FDind.LO %>% 
    reduce(inner_join, by = "Epoch")
  
  #close loop
})

stopImplicitCluster()

# Merge into 1 dataframe
res_df_LO <- res.LO %>% 
  bind_rows()

# Remove sp_richn, which is the same as nb_sp - redundant
res_df_LO$sp_richn<-NULL

# Format res for FD metrics
Res_FDmetrics_LO<- res_df_LO %>% 
  select(Epoch:fspe)
save(Res_FDmetrics_LO, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/LO_variation_metrics.RData")

# Melt data
FDmetrics_long_LO<- melt(Res_FDmetrics_LO, id.vars= "Epoch")
save(FDmetrics_long_LO, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/LO_variation_long_metrics.RData")

# Form dataframe of mean, median and standard deviation of all FD metrics
Taxon_LO <- Res_FDmetrics_LO %>% 
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

Taxon_LO$Epoch <- ordered(Taxon_CW$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                   "Pliocene","Pleistocene","Recent"))
Taxon_LO$Spp <- Spp_richness

# Save iteration data
save(Taxon_LO, file = "C:/Users/2022207/Dropbox/Jack's PhD/Chapter 2. FD changes over time/Analyses/Current Analyses/R codes/Taxon code/Data/Sensitivity tests/LO_Taxon_metrics.RData")



##################################################################################
# Plots
##################################################################################
# Plot empirical data against sensitivity tests
epochs_custom <- data.frame(
  name = c("Palaeocene","Eocene","Oligocene","Miocene","Pliocene","Pleistocene","Recent"), 
  max_age = c(1.5,2.5,3.5,4.5,5.5,6.5,7.6),
  min_age = c(0.4,1.5,2.5,3.5,4.5,5.5,6.5), 
  color = c("#FBA75F","#FDB46C","#FDC07A","#FFFF90","#FFFF99","#FFF2AE","#FEF2E0")
)

## Functional entities
FE_sens<-ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_line(data = Taxon_CH, aes(x= Epoch, y= FE_med, group = 1), color= "#F8766D", size = 1.5)+
  geom_line(data = Taxon_CW, aes(x= Epoch, y= FE_med, group = 1), color= "#B79F00", size = 1.5)+
  geom_line(data = Taxon_CE, aes(x= Epoch, y= FE_med, group = 1), color= "#00BA38", size = 1.5)+
  geom_line(data = Taxon_LC, aes(x= Epoch, y= FE_med, group = 1), color= "#00BFC4", size = 1.5)+
  geom_line(data = Taxon_XO, aes(x= Epoch, y= FE_med, group = 1), color= "#619CFF", size = 1.5)+
  geom_line(data = Taxon_XO, aes(x= Epoch, y= FE_med, group = 1), color= "#F564E3", size = 1.5)+
  geom_line(data=Taxon_var, aes(x=Epoch,y=FE_med, group=1), color="black", size=3)+
  labs(x = "", y = "# FEs")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional redundancy (with legend)
FRed_sens<-ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_line(data = Taxon_CH, aes(x= Epoch, y= Red_med, group = 1), color= "#F8766D", size = 1.5)+
  geom_line(data = Taxon_CW, aes(x= Epoch, y= Red_med, group = 1), color= "#B79F00", size = 1.5)+
  geom_line(data = Taxon_CE, aes(x= Epoch, y= Red_med, group = 1), color= "#00BA38", size = 1.5)+
  geom_line(data = Taxon_LC, aes(x= Epoch, y= Red_med, group = 1), color= "#00BFC4", size = 1.5)+
  geom_line(data = Taxon_XO, aes(x= Epoch, y= Red_med, group = 1), color= "#619CFF", size = 1.5)+
  geom_line(data = Taxon_XO, aes(x= Epoch, y= Red_med, group = 1), color= "#F564E3", size = 1.5)+
  geom_line(data=Taxon_var, aes(x=Epoch,y=Red_med, group=1), color="black", size=3)+
  labs(x = "", y = "FRed")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank()) +
  annotate("segment", xend = 6.3, yend = 6.5, x = 6, y = 6.5, color = "#F8766D", size = 1.5) +
  annotate("text", x = 6.32, y = 6.5, label = "No CH", hjust = 0, size = 5, color = "#F8766D") +
  annotate("segment", xend = 6.3, yend = 6.2, x = 6, y = 6.2, color = "#B79F00", size = 1.5) +
  annotate("text", x = 6.32, y = 6.2, label = "No CW", hjust = 0, size = 5, color = "#B79F00") +
  annotate("segment", xend = 6.3, yend = 5.9, x = 6, y = 5.9, color = "#00BA38", size = 1.5) +
  annotate("text", x = 6.32, y = 5.9, label = "No CE", hjust = 0, size = 5, color = "#00BA38") +
  annotate("segment", xend = 6.3, yend = 5.6, x = 6, y = 5.6, color = "#00BFC4", size = 1.5) +
  annotate("text", x = 6.32, y = 5.6, label = "No LC", hjust = 0, size = 5, color = "#00BFC4") +
  annotate("segment", xend = 6.3, yend = 5.3, x = 6, y = 5.3, color = "#619CFF", size = 1.5) +
  annotate("text", x = 6.32, y = 5.3, label = "No XO", hjust = 0, size = 5, color = "#619CFF") +
  annotate("segment", xend = 6.3, yend = 5, x = 6, y = 5, color = "#F564E3", size = 1.5) +
  annotate("text", x = 6.32, y = 5, label = "No LO", hjust = 0, size = 5, color = "#F564E3")+
  annotate("segment", xend = 6.3, yend = 4.7, x = 6, y = 4.7, color = "black", size = 1.5) +
  annotate("text", x = 6.32, y = 4.7, label = "Empirical", hjust = 0, size = 5, color = "black")

## Functional over-redundancy
FOred_sens<-ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_line(data = Taxon_CH, aes(x= Epoch, y= Ored_med, group = 1), color= "#F8766D", size = 1.5)+
  geom_line(data = Taxon_CW, aes(x= Epoch, y= Ored_med, group = 1), color= "#B79F00", size = 1.5)+
  geom_line(data = Taxon_CE, aes(x= Epoch, y= Ored_med, group = 1), color= "#00BA38", size = 1.5)+
  geom_line(data = Taxon_LC, aes(x= Epoch, y= Ored_med, group = 1), color= "#00BFC4", size = 1.5)+
  geom_line(data = Taxon_XO, aes(x= Epoch, y= Ored_med, group = 1), color= "#619CFF", size = 1.5)+
  geom_line(data = Taxon_XO, aes(x= Epoch, y= Ored_med, group = 1), color= "#F564E3", size = 1.5)+
  geom_line(data=Taxon_var, aes(x=Epoch,y=Ored_med, group=1), color="black", size=3)+
  labs(x = "", y = "FOred")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional richness
FRic_sens<-ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_line(data = Taxon_CH, aes(x= Epoch, y= FRic_med, group = 1), color= "#F8766D", size = 1.5)+
  geom_line(data = Taxon_CW, aes(x= Epoch, y= FRic_med, group = 1), color= "#B79F00", size = 1.5)+
  geom_line(data = Taxon_CE, aes(x= Epoch, y= FRic_med, group = 1), color= "#00BA38", size = 1.5)+
  geom_line(data = Taxon_LC, aes(x= Epoch, y= FRic_med, group = 1), color= "#00BFC4", size = 1.5)+
  geom_line(data = Taxon_XO, aes(x= Epoch, y= FRic_med, group = 1), color= "#619CFF", size = 1.5)+
  geom_line(data = Taxon_XO, aes(x= Epoch, y= FRic_med, group = 1), color= "#F564E3", size = 1.5)+
  geom_line(data=Taxon_var, aes(x=Epoch,y=FRic_med, group=1), color="black", size=3)+
  labs(x = "", y = "FRic")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional originality
FOri_sens<-ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_line(data = Taxon_CH, aes(x= Epoch, y= Fori_med, group = 1), color= "#F8766D", size = 1.5)+
  geom_line(data = Taxon_CW, aes(x= Epoch, y= Fori_med, group = 1), color= "#B79F00", size = 1.5)+
  geom_line(data = Taxon_CE, aes(x= Epoch, y= Fori_med, group = 1), color= "#00BA38", size = 1.5)+
  geom_line(data = Taxon_LC, aes(x= Epoch, y= Fori_med, group = 1), color= "#00BFC4", size = 1.5)+
  geom_line(data = Taxon_XO, aes(x= Epoch, y= Fori_med, group = 1), color= "#619CFF", size = 1.5)+
  geom_line(data = Taxon_XO, aes(x= Epoch, y= Fori_med, group = 1), color= "#F564E3", size = 1.5)+
  geom_line(data=Taxon_var, aes(x=Epoch,y=Fori_med, group=1), color="black", size=3)+
  labs(x = "", y = "FOri")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text.x = element_blank())

## Functional specialisation
FSpe_sens<-ggplot()+
  scale_x_discrete(limits=c("Palaeocene","Eocene","Oligocene","Miocene", "Pliocene", "Pleistocene","Recent"))+
  geom_line(data = Taxon_CH, aes(x= Epoch, y= Fspe_med, group = 1), color= "#F8766D", size = 1.5)+
  geom_line(data = Taxon_CW, aes(x= Epoch, y= Fspe_med, group = 1), color= "#B79F00", size = 1.5)+
  geom_line(data = Taxon_CE, aes(x= Epoch, y= Fspe_med, group = 1), color= "#00BA38", size = 1.5)+
  geom_line(data = Taxon_LC, aes(x= Epoch, y= Fspe_med, group = 1), color= "#00BFC4", size = 1.5)+
  geom_line(data = Taxon_XO, aes(x= Epoch, y= Fspe_med, group = 1), color= "#619CFF", size = 1.5)+
  geom_line(data = Taxon_XO, aes(x= Epoch, y= Fspe_med, group = 1), color= "#F564E3", size = 1.5)+
  geom_line(data=Taxon_var, aes(x=Epoch,y=Fspe_med, group=1), color="black", size=3)+
  labs(x = "Epoch", y = "FSpe")+
  theme_minimal() +
  theme_bw() +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size= 10, color= "black"),
        axis.title= element_text(size= 12), panel.background= element_rect(fill= "white"))+
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank())+
  coord_geo(
    dat = epochs_custom, pos = "bottom", expand = TRUE,
    skip = NULL, abbrv = TRUE, dat_is_discrete = TRUE,
    size = 5
  )

# Plot all metrics - produces Figure S9
plot_grid(FE_sens,FRed_sens,FOred_sens,FRic_sens,FOri_sens,FSpe_sens,
          labels= c("(a)","(b)","(c)","(d)","(e)","(f)"), 
          label_size = 12,align = "hv", label_fontface = "bold",  nrow=6)

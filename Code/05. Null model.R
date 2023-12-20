###################################################################################################################
# 05. Null model 
## This R code provides the null model used to assess functional diversity changes based on species richness 
## it produces the Rdata files needed to load null results for comparison to empirical analyses in Figure 2
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
library(xfun)

# Load data
load(file="~/Cleaned data.RData")

# Occurrence matrix
# Form FTU unique combinations, with epoch included
abun <- data %>% 
  group_by(Taxon_corrected) %>% 
  distinct(CH,CW,CE,LC,XO,LO,Epoch,Epoch_earliest,Epoch_latest,Current_status)

# Order taxa alphabetically like before and order Epochs
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

# Null model analyses
registerDoParallel(cores = 10)
getDoParWorkers()

n<-1000
res.Taxonvar.null<-NULL
res.Taxonvar.null<- lapply(1:n, function(x){
  sharks_traits.nullvar <- FTU_filtered %>%
    group_by(Taxon_corrected) %>%
    slice_sample(n=1, replace = FALSE)
  
  # Form species-trait matrix
  sharks_traits.nullvar <- sharks_traits.nullvar %>% 
    ungroup() %>%
    mutate(ID = c(1:537))
  sharks_traits.nullvar$ID <- numbers_to_words(sharks_traits.nullvar$ID)
  
  sharks_traits_null <- sharks_traits.nullvar %>% 
    column_to_rownames(var = "ID") %>% 
    select(-c(Taxon_corrected))
  
  
  # Form randomised number column in baskets
  baskets.null <- filtered_baskets.range %>%
    ungroup() %>% 
    mutate(vector = sample(1:537, 537, replace= FALSE))
  
  # Order 'vector' and set to column names
  baskets.null <- baskets.null[order(baskets.null$vector),]
  baskets.null$vector <- numbers_to_words(baskets.null$vector)
  
  baskets.null <- baskets.null %>%
    column_to_rownames(var = "vector") %>% 
    select(-c(Taxon_corrected,Current_status))
  
  # Transpose to form final occurrence matrix
  baskets_sharks_null <- t(baskets.null)
  
  # Make matrix and numeric
  baskets_sharks_weights_null <- data.matrix(baskets_sharks_null, rownames.force = NA)
  class(baskets_sharks_weights_null)<-"numeric"
  
  # Summarise & extract occurrences
  asb_sp_sharks_summ_null <- mFD::asb.sp.summary(asb_sp_w = baskets_sharks_weights_null)
  
  asb_sp_sharks_occ_null <- asb_sp_sharks_summ_null$asb_sp_occ
  
  ############################################################################################
  # Make trait_cat matrix
  sharks_traits_cat_null <- tibble(trait_name = c("CH", "CW", "CE", "LC", "XO", "LO"),
                                   trait_type = c("O", "O", "N", "N", "N", "N"),
                                   trait_weight = c(0.5,0.5,1,1,0.33,0.67))
  
  ############################################################################################
  # Form the functional entities
  sp_to_fe_sharks_null <- mFD::sp.to.fe(
    sp_tr       = sharks_traits_null, 
    tr_cat      = sharks_traits_cat_null, 
    fe_nm_type  = "fe_rank", 
    check_input = TRUE)
  
  # Calculate and display metrics
  alpha_fd_fe_sharks_null <- alpha.fd.fe(asb_sp_sharks_occ_null,
                                         sp_to_fe_sharks_null,
                                         ind_nm = c("fred", "fored", "fvuln"),
                                         check_input = TRUE,
                                         details_returned = TRUE)
  
  # Form dataframe
  FEmetrics.null <- as.data.frame(alpha_fd_fe_sharks_null$asb_fdfe) %>% 
    tibble::rownames_to_column("Epoch")
  FEmetrics.null$Epoch <- ordered(FEmetrics.null$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                                 "Pliocene","Pleistocene","Recent"))
  
  # Construct trait distance matrix using functional entities (fe_tr)
  sp_dist_sharks_null <- mFD::funct.dist(
    sp_tr         = sharks_traits_null,
    tr_cat        = sharks_traits_cat_null,
    metric        = "gower",
    scale_euclid  = "scale_center",
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  
  # Assess quality
  fspaces_quality_sharks_null <- mFD::quality.fspaces(
    sp_dist             = sp_dist_sharks_null,
    maxdim_pcoa         = 10,
    deviation_weighting = c("absolute","squared"),
    fdist_scaling       = FALSE,
    fdendro             = "average")
  
  # Return coordinates of axis
  sp_faxes_coord_sharks_null <- fspaces_quality_sharks_null$"details_fspaces"$"sp_pc_coord"
  
  # Form alpha diversity indices
  alpha_fd_indices_sharks_null <- mFD::alpha.fd.multidim(
    sp_faxes_coord   = sp_faxes_coord_sharks_null[ , c("PC1", "PC2", "PC3", "PC4")],
    asb_sp_w         = baskets_sharks_weights_null,
    ind_vect         = c("fric", "fori", "fspe"),
    scaling          = TRUE,
    check_input      = TRUE,
    details_returned = TRUE)
  
  fd_ind_values_sharks_null <- alpha_fd_indices_sharks_null$"functional_diversity_indices"
  
  ## Form dataframe
  FD.null <- as.data.frame(fd_ind_values_sharks_null) %>% 
    tibble::rownames_to_column("Epoch")
  FD.null$Epoch <- ordered(FD.null$Epoch, levels=c("Palaeocene","Eocene","Oligocene","Miocene",
                                                   "Pliocene","Pleistocene","Recent"))
  
  ###################################################################################
  # Form list to merge datasets
  FDind.null = list(FEmetrics.null,FD.null)
  
  # Output
  FDind.null %>% 
    reduce(inner_join, by = "Epoch")
  
  #close loop
})

# Merge lists into 1 dataframe
res_Taxonvar_null_list <- list(res.Taxonvar.null)

res_Taxonvar_null_df <- res_Taxonvar_null_list %>% 
  bind_rows()

# Remove sp_richn, which is the same as nb_sp - redundant
res_Taxonvar_null_df$sp_richn<-NULL

# Format res for FD metrics
Null_FDmetrics_taxonvar<- res_Taxonvar_null_df %>% 
  select(Epoch:fspe)
save(Null_FDmetrics_taxonvar, file = "~/Taxon_variation_null.RData")

# Melt data & save
FDmetrics_null_long_taxonvar<- melt(Null_FDmetrics_taxonvar, id.vars= "Epoch")
save(FDmetrics_null_long_taxonvar, file = "~/Taxon_variation_null_long.RData")

# Isolate individual functional diversity metrics
FEmetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "nb_fe")
FRedmetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fred")
FOredmetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fored")
FVulnmetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fvuln")
FRicmetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fric")
FOrimetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fori")
FSpemetrics_null_Taxonvar <- FDmetrics_null_long_taxonvar %>% 
  filter(variable == "fspe") 

# Test null model data for normal distribution - all show non-normal distribution (P < 0.05)
## Code below focuses on FE, FRed & FRic; all metrics per epoch due to high sample size for shapiro wilk test
FE_null_pal <- FEmetrics_null_Taxonvar %>% 
  filter(Epoch == "Palaeocene")
shapiro.test(FE_null_pal$value)
FRed_null_pal <- FRedmetrics_null_Taxonvar %>% 
  filter(Epoch == "Palaeocene")
shapiro.test(FRed_null_pal$value)
FRic_null_pal <- FRicmetrics_null_Taxonvar %>% 
  filter(Epoch == "Palaeocene")
shapiro.test(FRic_null_pal$value)

FE_null_eo <- FEmetrics_null_Taxonvar %>% 
  filter(Epoch == "Eocene")
shapiro.test(FE_null_eo$value)
FRed_null_eo <- FRedmetrics_null_Taxonvar %>% 
  filter(Epoch == "Eocene")
shapiro.test(FRed_null_eo$value)
FRic_null_eo <- FRicmetrics_null_Taxonvar %>% 
  filter(Epoch == "Eocene")
shapiro.test(FRic_null_eo$value)

FE_null_oli <- FEmetrics_null_Taxonvar %>% 
  filter(Epoch == "Oligocene")
shapiro.test(FE_null_oli$value)
FRed_null_oli <- FRedmetrics_null_Taxonvar %>% 
  filter(Epoch == "Oligocene")
shapiro.test(FRed_null_oli$value)
FRic_null_oli <- FRicmetrics_null_Taxonvar %>% 
  filter(Epoch == "Oligocene")
shapiro.test(FRic_null_oli$value)

FE_null_mio <- FEmetrics_null_Taxonvar %>% 
  filter(Epoch == "Miocene")
shapiro.test(FE_null_mio$value)
FRed_null_mio <- FRedmetrics_null_Taxonvar %>% 
  filter(Epoch == "Miocene")
shapiro.test(FRed_null_mio$value)
FRic_null_mio <- FRicmetrics_null_Taxonvar %>% 
  filter(Epoch == "Miocene")
shapiro.test(FRic_null_mio$value)

FE_null_plio <- FEmetrics_null_Taxonvar %>% 
  filter(Epoch == "Pliocene")
shapiro.test(FE_null_plio$value)
FRed_null_plio <- FRedmetrics_null_Taxonvar %>% 
  filter(Epoch == "Pliocene")
shapiro.test(FRed_null_plio$value)
FRic_null_plio <- FRicmetrics_null_Taxonvar %>% 
  filter(Epoch == "Pliocene")
shapiro.test(FRic_null_plio$value)

FE_null_ple <- FEmetrics_null_Taxonvar %>% 
  filter(Epoch == "Pleistocene")
shapiro.test(FE_null_ple$value)
FRed_null_ple <- FRedmetrics_null_Taxonvar %>% 
  filter(Epoch == "Pleistocene")
shapiro.test(FRed_null_ple$value)
FRic_null_ple <- FRicmetrics_null_Taxonvar %>% 
  filter(Epoch == "Pleistocene")
shapiro.test(FRic_null_ple$value)

FE_null_rec <- FEmetrics_null_Taxonvar %>% 
  filter(Epoch == "Recent")
shapiro.test(FE_null_rec$value)
FRed_null_rec <- FRedmetrics_null_Taxonvar %>% 
  filter(Epoch == "Recent")
shapiro.test(FRed_null_rec$value)
FRic_null_rec <- FRicmetrics_null_Taxonvar %>% 
  filter(Epoch == "Recent")
shapiro.test(FRic_null_rec$value)

